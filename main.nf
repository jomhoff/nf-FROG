#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion = true

/*
 * nf-FROG: Feedback-driven Reference Optimization through Genomic remapping.
 *
 * Biological invariants:
 *   - Iterative mapping references remain SNP-only and are never masked.
 *   - Only passing homozygous ALT SNPs alter intermediate references.
 *   - Final consensus masking is the union of reference ambiguity, gVCF
 *     callability failures, rejected variants, and optional external masks.
 *   - Every rejected variant site is masked rather than reverting to REF.
 *   - Samples retire independently once their iterative reference converges.
 */

params.samplesheet = null
params.indivs_file = null
params.readsdir = null
params.reads_suffix1 = '_R1.fastq.gz'
params.reads_suffix2 = '_R2.fastq.gz'
params.chromos_file = null
params.ref_file = null
params.ploidy_file = null
params.outputdir = 'results_nf_FROG'

params.max_rounds = 6
params.minimum_rounds = 2
params.convergence_max_changes = 0

params.trim_reads = false
params.fastp_extra = ''
params.mark_duplicates = true
params.kraken2_db = null
params.kraken2_extra = ''

params.gvcf_mode = 'BP_RESOLUTION'
params.min_gq = 20
params.min_qual = 30
params.min_qd = 2.0
params.min_mq = 40.0
params.max_snp_fs = 60.0
params.max_indel_fs = 200.0
params.max_snp_sor = 3.0
params.max_indel_sor = 10.0
params.min_read_pos_rank_sum = -8.0
params.min_mq_rank_sum = -12.5
params.min_het_allele_fraction = 0.2
params.max_het_major_fraction = 0.8
params.min_hom_alt_fraction = 0.9
params.min_allele_depth = 3
params.final_het_mode = 'iupac'
params.final_include_indels = false

params.min_depth_multiplier = 0.5
params.max_depth_multiplier = 2.0
params.depth_mad_multiplier = 6.0
params.depth_upper_quantile = 0.99
params.absolute_min_depth = 5
params.minimum_positive_depth_sites = 1000
params.mapq_min = 20
params.baseq_min = 20

params.repeat_mask_bed = null
params.low_mappability_bed = null
params.exclusion_bed = null
params.allow_reference_subset = false

params.haplotypecaller_extra = ''
params.genotypegvcfs_extra = ''

def parseTsv(inputFile) {
    def lines = file(inputFile, checkIfExists: true)
        .readLines()
        .findAll { it.trim() && !it.startsWith('#') }
    if (!lines) error "No records found in ${inputFile}"
    def header = lines[0].split('\t', -1) as List
    lines.drop(1).collect { line ->
        def values = line.split('\t', -1) as List
        header.withIndex().collectEntries { key, index ->
            [(key): index < values.size() ? values[index].trim() : '']
        }
    }
}

def getReferenceLengths(reference) {
    def lengths = [:].withDefault { 0L }
    def currentContig = null
    reference.eachLine { line ->
        if (line.startsWith('>')) {
            currentContig = line.substring(1).tokenize()[0]
        } else if (currentContig != null) {
            lengths[currentContig] += line.trim().size()
        }
    }
    lengths
}

process SUBSET_REFERENCE {
    label 'small'
    tag 'target chromosomes'

    publishDir "${params.outputdir}/00.reference", mode: 'copy'

    input:
    path reference
    val chromosome_list

    output:
    path 'target_reference.fa'

    script:
    """
    samtools faidx ${reference}
    samtools faidx ${reference} ${chromosome_list.collect { "'${it}'" }.join(' ')} > target_reference.fa
    test -s target_reference.fa
    """
}

process PREPROCESS_READ_UNIT {
    label 'preprocess'
    tag "${individual} ${unit}"

    publishDir "${params.outputdir}/00.read_qc/${individual}", mode: 'copy'

    input:
    tuple val(individual), val(unit), path(read1), path(read2),
          val(library), val(karyotype), val(default_ploidy), val(pcr_free)

    output:
    tuple val(individual), val(unit),
          path("${unit}.R1.fastq.gz"), path("${unit}.R2.fastq.gz"),
          val(library), val(karyotype), val(default_ploidy), val(pcr_free),
          path("${unit}.fastp.json"), path("${unit}.fastp.html")

    script:
    def qcOnly = params.trim_reads.toString().toBoolean()
        ? ''
        : '--disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g'
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${unit}.R1.fastq.gz \
        --out2 ${unit}.R2.fastq.gz \
        --json ${unit}.fastp.json \
        --html ${unit}.fastp.html \
        --thread ${task.cpus} \
        ${qcOnly} \
        ${params.fastp_extra}
    """
}

process KRAKEN_SCREEN {
    label 'preprocess'
    tag "${individual} ${unit}"

    publishDir "${params.outputdir}/00.read_qc/${individual}", mode: 'copy'

    input:
    tuple val(individual), val(unit), path(read1), path(read2)
    path kraken_db

    output:
    path "${unit}.kraken2.report"

    script:
    """
    kraken2 \
        --db ${kraken_db} \
        --paired ${read1} ${read2} \
        --threads ${task.cpus} \
        --report ${unit}.kraken2.report \
        --output /dev/null \
        ${params.kraken2_extra}
    """
}

process INDEX_REFERENCE {
    label 'mapping'
    tag "${reference.simpleName}"

    publishDir "${params.outputdir}/00.reference_indexes", mode: 'copy',
        saveAs: { filename -> "${reference.simpleName}/${filename}" }

    input:
    tuple val(ref_key), path(reference)

    output:
    tuple val(ref_key),
          path(reference),
          path("${reference}.bwt.2bit.64"),
          path("${reference}.0123"),
          path("${reference}.pac"),
          path("${reference}.ann"),
          path("${reference}.amb"),
          path("${reference}.fai"),
          path("${reference.simpleName}.dict")

    script:
    """
    samtools faidx ${reference}
    gatk CreateSequenceDictionary -R ${reference} -O ${reference.simpleName}.dict
    bwa-mem2 index ${reference}
    """
}

process MAP_READ_UNIT {
    label 'mapping'
    tag "${individual} ${unit} round${round}"

    input:
    tuple val(individual), val(unit), path(read1), path(read2), val(library),
          path(reference), path(bwt), path(zero123), path(pac), path(ann),
          path(amb), path(fai), path(dict), val(round)

    output:
    tuple val(individual), val(round), val(unit),
          path("${unit}.bam"), path("${unit}.bam.bai")

    script:
    """
    bwa-mem2 mem \
        -t ${task.cpus} \
        -R '@RG\\tID:${unit}\\tSM:${individual}\\tLB:${library}\\tPL:ILLUMINA' \
        ${reference} ${read1} ${read2} |
      samtools sort -@ ${task.cpus} -m 2G -o ${unit}.bam
    samtools index -@ ${task.cpus} ${unit}.bam
    """
}

process MERGE_AND_MARK_DUPLICATES {
    label 'mapping'
    tag "${individual} round${round}"

    publishDir "${params.outputdir}/01.bams/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), path(unit_bams)

    output:
    tuple val(individual), val(round),
          path("${individual}.round${round}.bam"),
          path("${individual}.round${round}.bam.bai"),
          path("${individual}.round${round}.duplicate_metrics.txt")

    script:
    def mergeCommand = unit_bams.size() == 1
        ? "cp ${unit_bams[0]} merged.bam"
        : "samtools merge -@ ${task.cpus} -f merged.bam ${unit_bams.join(' ')}"
    def markCommand = params.mark_duplicates.toString().toBoolean()
        ? """gatk MarkDuplicates \
              -I merged.bam \
              -O ${individual}.round${round}.bam \
              -M ${individual}.round${round}.duplicate_metrics.txt \
              --ASSUME_SORT_ORDER coordinate"""
        : """cp merged.bam ${individual}.round${round}.bam
             printf 'DUPLICATES_NOT_MARKED\\n' > ${individual}.round${round}.duplicate_metrics.txt"""
    """
    ${mergeCommand}
    ${markCommand}
    samtools index -@ ${task.cpus} ${individual}.round${round}.bam
    """
}

process ALIGNMENT_QC {
    label 'small'
    tag "${individual} round${round}"

    publishDir "${params.outputdir}/07.qc/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), path(bam), path(bai), path(reference), path(fai), path(dict)

    output:
    tuple val(individual), val(round),
          path("${individual}.round${round}.flagstat.txt"),
          path("${individual}.round${round}.alignment_summary.txt"),
          path("${individual}.round${round}.insert_size_metrics.txt"),
          path("${individual}.round${round}.insert_size.pdf")

    script:
    """
    samtools flagstat ${bam} > ${individual}.round${round}.flagstat.txt
    gatk CollectAlignmentSummaryMetrics \
        -R ${reference} -I ${bam} \
        -O ${individual}.round${round}.alignment_summary.txt
    gatk CollectInsertSizeMetrics \
        -I ${bam} \
        -O ${individual}.round${round}.insert_size_metrics.txt \
        -H ${individual}.round${round}.insert_size.pdf \
        --M 0.5
    """
}

process CALCULATE_DEPTH_THRESHOLDS {
    label 'small'
    tag "${individual}"

    publishDir "${params.outputdir}/06.coverage", mode: 'copy'

    input:
    tuple val(individual), val(round), path(bam), path(bai), path(fai)

    output:
    tuple val(individual), path("${individual}.round${round}.depth_thresholds.tsv")

    script:
    """
    samtools depth -aa -Q ${params.mapq_min} -q ${params.baseq_min} ${bam} |
      python3 "${projectDir}/bin/calculate_depth_thresholds.py" \
        --fai ${fai} \
        --min-multiplier ${params.min_depth_multiplier} \
        --max-multiplier ${params.max_depth_multiplier} \
        --mad-multiplier ${params.depth_mad_multiplier} \
        --upper-quantile ${params.depth_upper_quantile} \
        --absolute-minimum ${params.absolute_min_depth} \
        --minimum-positive-sites ${params.minimum_positive_depth_sites} \
        --output ${individual}.round${round}.depth_thresholds.tsv
    """
}

process HAPLOTYPECALLER_GVCF {
    label 'gatk'
    tag "${individual} round${round} ${interval} ploidy${ploidy}"

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          val(shard_order), val(interval), val(ploidy), val(pcr_free),
          path(bam), path(bai), path(reference), path(fai), path(dict)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome), val(shard_order),
          path("${chrom_order}.${shard_order}.g.vcf.gz"),
          path("${chrom_order}.${shard_order}.g.vcf.gz.tbi")

    script:
    def pcrArgument = pcr_free ? '--pcr-indel-model NONE' : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -L '${interval}' \
        --sample-ploidy ${ploidy} \
        -ERC ${params.gvcf_mode.toString().toUpperCase()} \
        --native-pair-hmm-threads ${task.cpus} \
        ${pcrArgument} \
        ${params.haplotypecaller_extra} \
        -O ${chrom_order}.${shard_order}.g.vcf.gz
    """
}

process GATHER_CHROMOSOME_GVCF {
    label 'gatk'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/02.gvcfs/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome), path(gvcfs)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.g.vcf.gz"),
          path("${chrom_order}.${chromosome}.g.vcf.gz.tbi")

    script:
    def inputs = gvcfs.sort { it.getName() }.collect { "-I ${it}" }.join(' ')
    """
    gatk GatherVcfs \
        ${inputs} \
        -O ${chrom_order}.${chromosome}.g.vcf.gz
    """
}

process GENOTYPE_GVCF {
    label 'gatk'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/02.variants/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(gvcf), path(gvcf_tbi), path(reference), path(fai), path(dict)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.genotyped.vcf.gz"),
          path("${chrom_order}.${chromosome}.genotyped.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \
        -R ${reference} \
        -V ${gvcf} \
        -L '${chromosome}' \
        ${params.genotypegvcfs_extra} \
        -O ${chrom_order}.${chromosome}.genotyped.vcf.gz
    """
}

process MAKE_GVCF_MASK {
    label 'small'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/08.masks/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(gvcf), path(gvcf_tbi), path(thresholds), path(fai)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.gvcf_mask.bed"),
          path("${chrom_order}.${chromosome}.gvcf_mask_stats.tsv")

    script:
    """
    python3 "${projectDir}/bin/gvcf_callability_mask.py" \
        --gvcf ${gvcf} \
        --fai ${fai} \
        --chromosome '${chromosome}' \
        --thresholds ${thresholds} \
        --minimum-gq ${params.min_gq} \
        --output ${chrom_order}.${chromosome}.gvcf_mask.bed \
        --stats ${chrom_order}.${chromosome}.gvcf_mask_stats.tsv
    """
}

process MAKE_PLOIDY_MASK {
    label 'small'
    tag "${individual} ${chromosome} ploidy mask"

    publishDir "${params.outputdir}/08.masks/${individual}/ploidy", mode: 'copy'

    input:
    tuple val(individual), val(chrom_order), val(chromosome), val(intervals)

    output:
    tuple val(individual), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.ploidy_mask.bed")

    script:
    def rows = intervals.collect { interval ->
        "${chromosome}\\t${interval.start - 1}\\t${interval.end}\\tPLOIDY_ZERO"
    }.join('\\n')
    """
    printf '%b' '${rows}${rows ? '\\n' : ''}' > ${chrom_order}.${chromosome}.ploidy_mask.bed
    """
}

process PREPARE_CONSENSUS_VCFS {
    label 'small'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/02.variants/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(vcf), path(vcf_tbi), path(thresholds), path(ploidy_mask)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.iterative.vcf.gz"),
          path("${chrom_order}.${chromosome}.iterative.vcf.gz.tbi"),
          path("${chrom_order}.${chromosome}.final.vcf.gz"),
          path("${chrom_order}.${chromosome}.final.vcf.gz.tbi"),
          path("${chrom_order}.${chromosome}.final_rejects.bed"),
          path("${chrom_order}.${chromosome}.filter_stats.tsv")

    script:
    def commonArguments = """
        --vcf ${vcf}
        --chromosome '${chromosome}'
        --thresholds ${thresholds}
        --minimum-gq ${params.min_gq}
        --minimum-qual ${params.min_qual}
        --minimum-qd ${params.min_qd}
        --minimum-mq ${params.min_mq}
        --maximum-snp-fs ${params.max_snp_fs}
        --maximum-indel-fs ${params.max_indel_fs}
        --maximum-snp-sor ${params.max_snp_sor}
        --maximum-indel-sor ${params.max_indel_sor}
        --minimum-read-pos-rank-sum ${params.min_read_pos_rank_sum}
        --minimum-mq-rank-sum ${params.min_mq_rank_sum}
        --minimum-het-allele-fraction ${params.min_het_allele_fraction}
        --maximum-het-major-fraction ${params.max_het_major_fraction}
        --minimum-hom-alt-fraction ${params.min_hom_alt_fraction}
        --minimum-allele-depth ${params.min_allele_depth}
        --exclude-bed ${ploidy_mask}
    """.stripIndent().replaceAll('\\n', ' ')
    """
    python3 "${projectDir}/bin/prepare_consensus_vcf.py" \
        ${commonArguments} \
        --het-mode retain \
        --output-vcf iterative.vcf \
        --reject-bed iterative.rejects.bed \
        --stats iterative.stats.tsv
    bgzip -c iterative.vcf > ${chrom_order}.${chromosome}.iterative.vcf.gz
    tabix -f -p vcf ${chrom_order}.${chromosome}.iterative.vcf.gz

    python3 "${projectDir}/bin/prepare_consensus_vcf.py" \
        ${commonArguments} \
        --het-mode ${params.final_het_mode.toString().toLowerCase()} \
        --output-vcf final.vcf \
        --reject-bed ${chrom_order}.${chromosome}.final_rejects.bed \
        --stats ${chrom_order}.${chromosome}.filter_stats.tsv
    bgzip -c final.vcf > ${chrom_order}.${chromosome}.final.vcf.gz
    tabix -f -p vcf ${chrom_order}.${chromosome}.final.vcf.gz
    """
}

process VARIANT_QC {
    label 'small'
    tag "${individual} round${round} ${chromosome} variant QC"

    publishDir "${params.outputdir}/07.qc/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(vcf), path(vcf_tbi)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.variant_qc.tsv")

    script:
    """
    python3 "${projectDir}/bin/variant_qc_metrics.py" \
        --vcf ${vcf} \
        --chromosome '${chromosome}' \
        --output ${chrom_order}.${chromosome}.variant_qc.tsv
    """
}

process MAKE_REFERENCE_MASK {
    label 'small'
    tag "${individual} round${round} ${chromosome}"

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome), path(reference)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.reference_mask.bed")

    script:
    """
    python3 "${projectDir}/bin/fasta_n_mask.py" \
        --fasta ${reference} \
        --chromosome '${chromosome}' \
        --output ${chrom_order}.${chromosome}.reference_mask.bed
    """
}

process MERGE_FINAL_MASKS {
    label 'small'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/08.masks/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(gvcf_mask), path(variant_rejects), path(reference_mask),
          path(ploidy_mask),
          path(repeat_mask), path(mappability_mask), path(exclusion_mask)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.final_mask.bed"),
          path("${chrom_order}.${chromosome}.final_mask.reasons.bed"),
          path("${chrom_order}.${chromosome}.final_mask.reason_*.bed")

    script:
    """
    python3 "${projectDir}/bin/merge_reason_masks.py" \
        --chromosome '${chromosome}' \
        --mask GVCF=${gvcf_mask} \
        --mask VARIANT=${variant_rejects} \
        --mask REFERENCE=${reference_mask} \
        --mask PLOIDY=${ploidy_mask} \
        --mask REPEAT=${repeat_mask} \
        --mask LOW_MAPPABILITY=${mappability_mask} \
        --mask EXCLUSION=${exclusion_mask} \
        --prefix ${chrom_order}.${chromosome}.final_mask
    """
}

process BUILD_ITERATIVE_CHROMOSOME {
    label 'small'
    tag "${individual} round${round} ${chromosome}"

    publishDir "${params.outputdir}/05.round_consensus/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(vcf), path(vcf_tbi), path(reference), path(fai)

    output:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.iterative.fa")

    script:
    """
    samtools faidx ${reference} '${chromosome}' > chromosome.fa
    bcftools consensus \
        -f chromosome.fa \
        -s '${individual}' \
        -H I \
        ${vcf} \
        -o ${chrom_order}.${chromosome}.iterative.fa
    """
}

process BUILD_ITERATIVE_REFERENCE {
    label 'small'
    tag "${individual} round${round}"

    publishDir "${params.outputdir}/04.consensus_refs/${individual}", mode: 'copy'

    input:
    tuple val(individual), val(round), path(chromosome_fastas)

    output:
    tuple val(individual), val(round), path("${individual}.round${round}.iterative.fa")

    script:
    """
    for fasta in \$(printf '%s\\n' ${chromosome_fastas} | sort -V); do
        cat "\$fasta"
    done > ${individual}.round${round}.iterative.fa
    test -s ${individual}.round${round}.iterative.fa
    """
}

process COMPARE_REFERENCES {
    label 'small'
    tag "${individual} round${round}"

    publishDir "${params.outputdir}/07.qc/${individual}/round${round}", mode: 'copy'

    input:
    tuple val(individual), val(round),
          path(previous_reference), path(current_reference), path(new_reference)

    output:
    tuple val(individual), val(round),
          path("${individual}.round${round}.next.fa"),
          path("${individual}.round${round}.current.fa"),
          path("${individual}.round${round}.convergence.tsv"),
          path("${individual}.round${round}.status.txt")

    script:
    """
    cp ${new_reference} ${individual}.round${round}.next.fa
    cp ${current_reference} ${individual}.round${round}.current.fa
    python3 "${projectDir}/bin/compare_references.py" \
        --previous ${previous_reference} \
        --current ${current_reference} \
        --new ${new_reference} \
        --round ${round} \
        --minimum-rounds ${params.minimum_rounds} \
        --maximum-changes ${params.convergence_max_changes} \
        --metrics ${individual}.round${round}.convergence.tsv \
        --status ${individual}.round${round}.status.txt
    """
}

process BUILD_FINAL_CHROMOSOME {
    label 'small'
    tag "${individual} final ${chromosome}"

    publishDir "${params.outputdir}/03.consensus/${individual}/chromosomes", mode: 'copy'

    input:
    tuple val(individual), val(round), val(chrom_order), val(chromosome),
          path(vcf), path(vcf_tbi), path(mask), path(reference), path(fai)

    output:
    tuple val(individual), val(chrom_order), val(chromosome),
          path("${chrom_order}.${chromosome}.consensus.fa")

    script:
    """
    samtools faidx ${reference} '${chromosome}' > chromosome.unmasked.fa
    bedtools maskfasta \
        -fi chromosome.unmasked.fa \
        -bed ${mask} \
        -fo chromosome.masked.fa
    samtools faidx chromosome.masked.fa
    bcftools consensus \
        -f chromosome.masked.fa \
        -s '${individual}' \
        -H I \
        ${vcf} \
        -o ${chrom_order}.${chromosome}.consensus.fa
    """
}

process BUILD_FINAL_REFERENCE {
    label 'small'
    tag "${individual}"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(chromosome_fastas)

    output:
    tuple val(individual),
          path("${individual}.consensus.fa"),
          path("${individual}.consensus.fa.fai")

    script:
    """
    for fasta in \$(printf '%s\\n' ${chromosome_fastas} | sort -V); do
        cat "\$fasta"
    done > ${individual}.consensus.fa
    samtools faidx ${individual}.consensus.fa
    """
}

process GATHER_FINAL_VCF {
    label 'small'
    tag "${individual} final VCF"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(vcfs)

    output:
    tuple val(individual),
          path("${individual}.final.current_reference.vcf.gz"),
          path("${individual}.final.current_reference.vcf.gz.tbi")

    script:
    """
    bcftools concat -a -Oz \
        -o ${individual}.final.current_reference.vcf.gz \
        ${vcfs.sort { it.getName() }.join(' ')}
    tabix -f -p vcf ${individual}.final.current_reference.vcf.gz
    """
}

process GATHER_FINAL_GVCF {
    label 'gatk'
    tag "${individual} final gVCF"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(gvcfs)

    output:
    tuple val(individual),
          path("${individual}.final.g.vcf.gz"),
          path("${individual}.final.g.vcf.gz.tbi")

    script:
    def inputs = gvcfs.sort { it.getName() }.collect { "-I ${it}" }.join(' ')
    """
    gatk GatherVcfs ${inputs} -O ${individual}.final.g.vcf.gz
    """
}

process GATHER_FINAL_MASK {
    label 'small'
    tag "${individual} final mask"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(masks), path(reason_masks)

    output:
    tuple val(individual),
          path("${individual}.noncallable_mask.bed"),
          path("${individual}.noncallable_mask.reasons.bed")

    script:
    """
    for bed in \$(printf '%s\\n' ${masks} | sort -V); do
        cat "\$bed"
    done > ${individual}.noncallable_mask.bed
    for bed in \$(printf '%s\\n' ${reason_masks} | sort -V); do
        cat "\$bed"
    done > ${individual}.noncallable_mask.reasons.bed
    """
}

process MAKE_CALLABLE_BED {
    label 'small'
    tag "${individual} callable BED"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(consensus), path(fai), path(noncallable_mask)

    output:
    tuple val(individual), path("${individual}.callable.bed")

    script:
    """
    bedtools complement \
        -i ${noncallable_mask} \
        -g ${fai} \
        > ${individual}.callable.bed
    """
}

process ORIGINAL_REFERENCE_DIFF {
    label 'small'
    tag "${individual} original-coordinate variants"

    publishDir "${params.outputdir}/03.consensus/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(consensus), path(consensus_fai), path(original_reference)

    output:
    tuple val(individual),
          path("${individual}.original_reference.vcf.gz"),
          path("${individual}.original_reference.vcf.gz.tbi"),
          path("${individual}.original_reference_diff.tsv")

    script:
    """
    python3 "${projectDir}/bin/fasta_diff_to_vcf.py" \
        --reference ${original_reference} \
        --consensus ${consensus} \
        --sample '${individual}' \
        --output ${individual}.original_reference.vcf \
        --stats ${individual}.original_reference_diff.tsv
    bgzip -c ${individual}.original_reference.vcf > ${individual}.original_reference.vcf.gz
    tabix -f -p vcf ${individual}.original_reference.vcf.gz
    """
}

process FINAL_SAMPLE_QC {
    label 'small'
    tag "${individual} final QC"

    publishDir "${params.outputdir}/07.qc/${individual}", mode: 'copy'

    input:
    tuple val(individual), path(consensus), path(fai), path(mask), path(vcf)

    output:
    path "${individual}.final_qc.tsv"

    script:
    """
    total_bases=\$(awk '{sum+=\$2} END{print sum+0}' ${fai})
    masked_bases=\$(awk '{sum+=\$3-\$2} END{print sum+0}' ${mask})
    variants=\$(bcftools view -H ${vcf} | wc -l | tr -d ' ')
    callable_bases=\$((total_bases-masked_bases))
    awk -v sample='${individual}' \
        -v total="\$total_bases" \
        -v masked="\$masked_bases" \
        -v callable="\$callable_bases" \
        -v variants="\$variants" \
        'BEGIN {
            OFS="\\t";
            print "sample","total_bases","callable_bases","masked_bases","callable_fraction","final_variants";
            print sample,total,callable,masked,(total ? callable/total : 0),variants
        }' > ${individual}.final_qc.tsv
    """
}

process CONVERGENCE_SUMMARY {
    label 'small'
    tag 'all convergence rounds'

    publishDir "${params.outputdir}/07.qc", mode: 'copy'

    input:
    path metrics

    output:
    path 'all_samples.convergence.tsv'

    script:
    """
    printf 'sample\\tround\\tchanged_bases\\tchanged_acgt_bases\\toscillating_bases\\tconverged\\n' > all_samples.convergence.tsv
    for file in ${metrics}; do
        sample=\$(basename "\$file" | sed 's/\\.round[0-9]*\\.convergence\\.tsv//')
        tail -n +2 "\$file" | awk -v sample="\$sample" 'BEGIN{OFS="\\t"} {print sample,\$0}'
    done | sort -k1,1 -k2,2n >> all_samples.convergence.tsv
    """
}

workflow CALL_ROUND {
    take:
    processed_reads_ch
    active_refs_ch
    targets_ch
    sample_meta_ch
    chromosome_meta_ch
    ploidy_masks_ch
    round
    max_rounds
    repeat_mask_ch
    mappability_mask_ch
    exclusion_mask_ch

    main:
    refs_to_index_ch = active_refs_ch
        .map { ref_key, individual, current_ref, previous_ref -> tuple(ref_key, current_ref) }
        .distinct()
    INDEX_REFERENCE(refs_to_index_ch)

    indexed_refs_ch = active_refs_ch
        .combine(INDEX_REFERENCE.out, by: 0)
        .map { ref_key, individual, current_ref, previous_ref,
               indexed_ref, bwt, zero123, pac, ann, amb, fai, dict ->
            tuple(individual, current_ref, previous_ref, indexed_ref, bwt, zero123, pac, ann, amb, fai, dict)
        }

    map_input_ch = processed_reads_ch
        .combine(indexed_refs_ch, by: 0)
        .map { individual, unit, read1, read2, library, karyotype, default_ploidy, pcr_free,
               fastp_json, fastp_html, current_ref, previous_ref, indexed_ref,
               bwt, zero123, pac, ann, amb, fai, dict ->
            tuple(individual, unit, read1, read2, library,
                  indexed_ref, bwt, zero123, pac, ann, amb, fai, dict, round)
        }
    MAP_READ_UNIT(map_input_ch)

    merged_bam_input_ch = MAP_READ_UNIT.out
        .groupTuple(by: [0, 1])
        .map { individual, rnd, units, bams, bais ->
            tuple(individual, rnd, bams.flatten())
        }
    MERGE_AND_MARK_DUPLICATES(merged_bam_input_ch)

    ref_tools_ch = indexed_refs_ch
        .map { individual, current_ref, previous_ref, indexed_ref, bwt, zero123, pac, ann, amb, fai, dict ->
            tuple(individual, indexed_ref, fai, dict, current_ref, previous_ref)
        }

    bam_ref_ch = MERGE_AND_MARK_DUPLICATES.out
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, bam, bai, duplicate_metrics, reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, bam, bai, reference, fai, dict)
        }
    ALIGNMENT_QC(bam_ref_ch)

    threshold_input_ch = MERGE_AND_MARK_DUPLICATES.out
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, bam, bai, duplicate_metrics, reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, bam, bai, fai)
        }
    active_thresholds_ch = CALCULATE_DEPTH_THRESHOLDS(threshold_input_ch)

    call_input_ch = MERGE_AND_MARK_DUPLICATES.out
        .combine(targets_ch, by: 0)
        .combine(sample_meta_ch, by: 0)
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, bam, bai, duplicate_metrics,
               chrom_order, chromosome, shard_order, interval, ploidy,
               karyotype, default_ploidy, pcr_free,
               reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, chrom_order, chromosome, shard_order, interval,
                  ploidy, pcr_free, bam, bai, reference, fai, dict)
        }
    HAPLOTYPECALLER_GVCF(call_input_ch)

    gathered_gvcf_input_ch = HAPLOTYPECALLER_GVCF.out
        .groupTuple(by: [0, 1, 2, 3])
        .map { individual, rnd, chrom_order, chromosome, shard_orders, gvcfs, tbis ->
            tuple(individual, rnd, chrom_order, chromosome, gvcfs.flatten())
        }
    GATHER_CHROMOSOME_GVCF(gathered_gvcf_input_ch)

    genotype_input_ch = GATHER_CHROMOSOME_GVCF.out
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, chrom_order, chromosome, gvcf, gvcf_tbi,
               reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, chrom_order, chromosome, gvcf, gvcf_tbi, reference, fai, dict)
        }
    GENOTYPE_GVCF(genotype_input_ch)
    VARIANT_QC(GENOTYPE_GVCF.out)

    filter_input_ch = GENOTYPE_GVCF.out
        .combine(active_thresholds_ch, by: 0)
        .combine(ploidy_masks_ch.map { individual, chrom_order, chromosome, mask ->
            tuple(individual, round, chrom_order, chromosome, mask)
        }, by: [0, 1, 2, 3])
        .map { individual, rnd, chrom_order, chromosome, vcf, vcf_tbi, thresholds, ploidy_mask ->
            tuple(individual, rnd, chrom_order, chromosome, vcf, vcf_tbi, thresholds, ploidy_mask)
        }
    PREPARE_CONSENSUS_VCFS(filter_input_ch)

    gvcf_mask_input_ch = GATHER_CHROMOSOME_GVCF.out
        .combine(active_thresholds_ch, by: 0)
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, chrom_order, chromosome, gvcf, gvcf_tbi,
               thresholds, reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, chrom_order, chromosome, gvcf, gvcf_tbi, thresholds, fai)
        }
    MAKE_GVCF_MASK(gvcf_mask_input_ch)

    reference_mask_input_ch = chromosome_meta_ch
        .combine(ref_tools_ch)
        .map { chrom_order, chromosome, individual, reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, round, chrom_order, chromosome, reference)
        }
    MAKE_REFERENCE_MASK(reference_mask_input_ch)

    final_mask_input_ch = MAKE_GVCF_MASK.out
        .combine(PREPARE_CONSENSUS_VCFS.out, by: [0, 1, 2, 3])
        .combine(MAKE_REFERENCE_MASK.out, by: [0, 1, 2, 3])
        .combine(ploidy_masks_ch.map { individual, chrom_order, chromosome, mask ->
            tuple(individual, round, chrom_order, chromosome, mask)
        }, by: [0, 1, 2, 3])
        .combine(repeat_mask_ch)
        .combine(mappability_mask_ch)
        .combine(exclusion_mask_ch)
        .map { individual, rnd, chrom_order, chromosome, gvcf_mask, gvcf_stats,
               iterative_vcf, iterative_tbi, final_vcf, final_tbi, variant_rejects, filter_stats,
               reference_mask, ploidy_mask, repeat_mask, mappability_mask, exclusion_mask ->
            tuple(individual, rnd, chrom_order, chromosome, gvcf_mask, variant_rejects,
                  reference_mask, ploidy_mask, repeat_mask, mappability_mask, exclusion_mask)
        }
    MERGE_FINAL_MASKS(final_mask_input_ch)

    iterative_consensus_input_ch = PREPARE_CONSENSUS_VCFS.out
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, chrom_order, chromosome,
               iterative_vcf, iterative_tbi, final_vcf, final_tbi, variant_rejects, filter_stats,
               reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, chrom_order, chromosome, iterative_vcf, iterative_tbi, reference, fai)
        }
    BUILD_ITERATIVE_CHROMOSOME(iterative_consensus_input_ch)

    new_refs_ch = BUILD_ITERATIVE_CHROMOSOME.out
        .groupTuple(by: [0, 1])
        .map { individual, rnd, chrom_orders, chromosomes, fastas ->
            tuple(individual, rnd, fastas.flatten())
        }
        | BUILD_ITERATIVE_REFERENCE

    compare_input_ch = active_refs_ch
        .map { ref_key, individual, current_ref, previous_ref -> tuple(individual, previous_ref, current_ref) }
        .combine(new_refs_ch.map { individual, rnd, new_ref -> tuple(individual, rnd, new_ref) }, by: 0)
        .map { individual, previous_ref, current_ref, rnd, new_ref ->
            tuple(individual, rnd, previous_ref, current_ref, new_ref)
        }
    COMPARE_REFERENCES(compare_input_ch)

    status_ch = COMPARE_REFERENCES.out
        .map { individual, rnd, next_ref, current_ref, metrics, status_file ->
            tuple(individual, rnd, status_file.text.trim(), next_ref, current_ref, metrics)
        }

    next_active_refs_ch = status_ch
        .filter { individual, rnd, status, next_ref, current_ref, metrics ->
            status == 'active' && rnd < max_rounds
        }
        .map { individual, rnd, status, next_ref, current_ref, metrics ->
            tuple(next_ref.getName(), individual, next_ref, current_ref)
        }

    terminal_refs_ch = status_ch
        .filter { individual, rnd, status, next_ref, current_ref, metrics ->
            status == 'converged' || rnd >= max_rounds
        }
        .map { individual, rnd, status, next_ref, current_ref, metrics ->
            tuple(next_ref.getName(), individual, next_ref)
        }

    terminal_keys_ch = status_ch
        .filter { individual, rnd, status, next_ref, current_ref, metrics ->
            status == 'converged' || rnd >= max_rounds
        }
        .map { individual, rnd, status, next_ref, current_ref, metrics -> tuple(individual, rnd) }

    round_artifacts_ch = PREPARE_CONSENSUS_VCFS.out
        .combine(MERGE_FINAL_MASKS.out, by: [0, 1, 2, 3])
        .combine(GATHER_CHROMOSOME_GVCF.out, by: [0, 1, 2, 3])
        .combine(GENOTYPE_GVCF.out, by: [0, 1, 2, 3])
        .combine(ref_tools_ch, by: 0)
        .map { individual, rnd, chrom_order, chromosome,
               iterative_vcf, iterative_tbi, final_vcf, final_tbi, variant_rejects, filter_stats,
               final_mask, reason_mask, reason_files,
               gvcf, gvcf_tbi, genotyped_vcf, genotyped_tbi,
               reference, fai, dict, current_ref, previous_ref ->
            tuple(individual, rnd, chrom_order, chromosome,
                  final_vcf, final_tbi, final_mask, reason_mask,
                  gvcf, gvcf_tbi, genotyped_vcf, genotyped_tbi,
                  reference, fai, dict, filter_stats)
        }

    terminal_artifacts_ch = round_artifacts_ch
        .combine(terminal_keys_ch, by: [0, 1])

    emit:
    next_active_refs = next_active_refs_ch
    terminal_refs = terminal_refs_ch
    terminal_artifacts = terminal_artifacts_ch
    convergence_metrics = COMPARE_REFERENCES.out.map { it[4] }
}

workflow RUN_ITERATIONS {
    take:
    processed_reads_ch
    active_refs_ch
    targets_ch
    sample_meta_ch
    chromosome_meta_ch
    ploidy_masks_ch
    current_round
    max_rounds
    terminal_refs_acc
    terminal_artifacts_acc
    convergence_metrics_acc
    repeat_mask_ch
    mappability_mask_ch
    exclusion_mask_ch

    main:
    CALL_ROUND(
        processed_reads_ch,
        active_refs_ch,
        targets_ch,
        sample_meta_ch,
        chromosome_meta_ch,
        ploidy_masks_ch,
        current_round,
        max_rounds,
        repeat_mask_ch,
        mappability_mask_ch,
        exclusion_mask_ch
    )

    def accumulatedTerminalRefs = terminal_refs_acc.mix(CALL_ROUND.out.terminal_refs)
    def accumulatedArtifacts = terminal_artifacts_acc.mix(CALL_ROUND.out.terminal_artifacts)
    def accumulatedMetrics = convergence_metrics_acc.mix(CALL_ROUND.out.convergence_metrics)

    if (current_round < max_rounds) {
        RUN_ITERATIONS(
            processed_reads_ch,
            CALL_ROUND.out.next_active_refs,
            targets_ch,
            sample_meta_ch,
            chromosome_meta_ch,
            ploidy_masks_ch,
            current_round + 1,
            max_rounds,
            accumulatedTerminalRefs,
            accumulatedArtifacts,
            accumulatedMetrics,
            repeat_mask_ch,
            mappability_mask_ch,
            exclusion_mask_ch
        )
        finalRefs = RUN_ITERATIONS.out.final_refs
        finalArtifacts = RUN_ITERATIONS.out.final_artifacts
        allMetrics = RUN_ITERATIONS.out.convergence_metrics
    } else {
        finalRefs = accumulatedTerminalRefs
        finalArtifacts = accumulatedArtifacts
        allMetrics = accumulatedMetrics
    }

    emit:
    final_refs = finalRefs
    final_artifacts = finalArtifacts
    convergence_metrics = allMetrics
}

workflow {
    if (!params.ref_file || !params.chromos_file) {
        error "Required parameters: --ref_file and --chromos_file"
    }
    if (!params.samplesheet && (!params.indivs_file || !params.readsdir)) {
        error "Provide --samplesheet, or legacy --indivs_file plus --readsdir"
    }
    if (params.final_include_indels.toString().toBoolean()) {
        error "Indels are intentionally disabled so FASTAs retain original coordinates; use the delivered VCF for indels."
    }

    maxRounds = params.max_rounds as int
    minimumRounds = params.minimum_rounds as int
    if (maxRounds < 1 || minimumRounds < 1 || minimumRounds > maxRounds) {
        error "Require 1 <= --minimum_rounds <= --max_rounds"
    }
    finalHetMode = params.final_het_mode.toString().toLowerCase()
    if (!(finalHetMode in ['iupac', 'major', 'mask'])) {
        error "--final_het_mode must be iupac, major, or mask"
    }
    gvcfMode = params.gvcf_mode.toString().toUpperCase()
    if (!(gvcfMode in ['GVCF', 'BP_RESOLUTION'])) {
        error "--gvcf_mode must be GVCF or BP_RESOLUTION"
    }

    referenceFile = file(params.ref_file, checkIfExists: true)
    targetChromosomes = file(params.chromos_file, checkIfExists: true)
        .readLines()
        .collect { it.trim() }
        .findAll { it && !it.startsWith('#') }
    if (!targetChromosomes) {
        error "No chromosomes found in ${params.chromos_file}"
    }

    referenceLengths = getReferenceLengths(referenceFile)
    absentChromosomes = targetChromosomes.findAll { !referenceLengths.containsKey(it) }
    if (absentChromosomes) {
        error "Chromosomes absent from reference: ${absentChromosomes.join(', ')}"
    }
    omittedChromosomes = referenceLengths.keySet().findAll { !targetChromosomes.contains(it) }
    if (omittedChromosomes && !params.allow_reference_subset.toString().toBoolean()) {
        error "chromos_file omits ${omittedChromosomes.size()} reference contigs. Use --allow_reference_subset true only for an intentional subset."
    }

    readUnits = []
    if (params.samplesheet) {
        rows = parseTsv(params.samplesheet)
        ['individual', 'read1', 'read2'].each { column ->
            if (!rows[0].containsKey(column)) {
                error "samplesheet is missing '${column}'"
            }
        }
        rows.eachWithIndex { row, index ->
            def sampleId = row.individual
            def unitId = row.read_group ?: "${sampleId}.unit${index + 1}"
            readUnits << [
                individual: sampleId,
                unit: unitId,
                read1: file(row.read1, checkIfExists: true),
                read2: file(row.read2, checkIfExists: true),
                library: row.library ?: sampleId,
                karyotype: row.karyotype ?: 'default',
                default_ploidy: (row.default_ploidy ?: '2') as int,
                pcr_free: (row.pcr_free ?: 'false').toBoolean()
            ]
        }
    } else {
        file(params.indivs_file, checkIfExists: true)
            .readLines()
            .collect { it.trim() }
            .findAll { it && !it.startsWith('#') }
            .eachWithIndex { legacySampleId, index ->
                readUnits << [
                    individual: legacySampleId,
                    unit: "${legacySampleId}.unit${index + 1}",
                    read1: file("${params.readsdir}/${legacySampleId}${params.reads_suffix1}", checkIfExists: true),
                    read2: file("${params.readsdir}/${legacySampleId}${params.reads_suffix2}", checkIfExists: true),
                    library: legacySampleId,
                    karyotype: 'default',
                    default_ploidy: 2,
                    pcr_free: false
                ]
            }
    }

    samples = readUnits.groupBy { it.individual }.collect { sampleId, units ->
        def karyotypes = units.collect { it.karyotype }.unique()
        def ploidies = units.collect { it.default_ploidy }.unique()
        def pcrStates = units.collect { it.pcr_free }.unique()
        if (karyotypes.size() != 1 || ploidies.size() != 1 || pcrStates.size() != 1) {
            error "Conflicting sample metadata across read units for ${sampleId}"
        }
        [
            individual: sampleId,
            karyotype: karyotypes[0],
            default_ploidy: ploidies[0],
            pcr_free: pcrStates[0]
        ]
    }
    if (!samples) {
        error "No samples were defined"
    }
    safeIdentifier = ~/^[A-Za-z0-9_.-]+$/
    readUnits.each { row ->
        if (!(row.individual ==~ safeIdentifier) ||
            !(row.unit ==~ safeIdentifier) ||
            !(row.library ==~ safeIdentifier)) {
            error "individual, read_group, and library values may contain only letters, numbers, '.', '_', and '-'"
        }
    }
    if (readUnits.collect { it.unit }.size() != readUnits.collect { it.unit }.unique().size()) {
        error "read_group values must be unique across the sample sheet"
    }
    targetChromosomes.each { targetChromosome ->
        if (!(targetChromosome ==~ safeIdentifier)) {
            error "Chromosome names may contain only letters, numbers, '.', '_', and '-'"
        }
    }

    ploidyRows = params.ploidy_file ? parseTsv(params.ploidy_file) : []
    if (ploidyRows) {
        ['karyotype', 'chromosome', 'start', 'end', 'ploidy'].each { column ->
            if (!ploidyRows[0].containsKey(column)) {
                error "ploidy_file is missing '${column}'"
            }
        }
    }

    callTargets = []
    samples.each { sample ->
        targetChromosomes.eachWithIndex { targetChromosome, chromosomeIndex ->
            def rules = ploidyRows
                .findAll { it.karyotype == sample.karyotype && it.chromosome == targetChromosome }
                .collect {
                    [
                        start: it.start as int,
                        end: it.end == '*' ? referenceLengths[targetChromosome] as int : it.end as int,
                        ploidy: it.ploidy as int
                    ]
                }
                .sort { it.start }
            if (!rules) {
                rules = [[start: 1, end: referenceLengths[targetChromosome] as int, ploidy: sample.default_ploidy]]
            }
            def expectedStart = 1
            rules.eachWithIndex { rule, shardIndex ->
                if (rule.ploidy < 0 || rule.start != expectedStart || rule.end < rule.start || rule.end > referenceLengths[targetChromosome]) {
                    error "Ploidy intervals must form a complete, non-overlapping partition: ${sample.individual} ${targetChromosome}"
                }
                callTargets << [
                    individual: sample.individual,
                    chrom_order: String.format('%06d', chromosomeIndex),
                    chromosome: targetChromosome,
                    shard_order: String.format('%06d', shardIndex),
                    interval: "${targetChromosome}:${rule.start}-${rule.end}",
                    start: rule.start,
                    end: rule.end,
                    ploidy: Math.max(1, rule.ploidy),
                    declared_ploidy: rule.ploidy
                ]
                expectedStart = rule.end + 1
            }
            if (expectedStart != referenceLengths[targetChromosome] + 1) {
                error "Ploidy intervals do not cover all of ${targetChromosome} for ${sample.individual}"
            }
        }
    }

    ploidyMaskDefs = callTargets
        .groupBy { "${it.individual}\t${it.chrom_order}\t${it.chromosome}" }
        .collect { key, targets ->
            def parts = key.split('\t')
            [
                individual: parts[0],
                chrom_order: parts[1],
                chromosome: parts[2],
                intervals: targets
                    .findAll { it.declared_ploidy == 0 }
                    .collect { [start: it.start, end: it.end] }
            ]
        }

    repeatMask = file(params.repeat_mask_bed ?: "${projectDir}/assets/empty.bed", checkIfExists: true)
    mappabilityMask = file(params.low_mappability_bed ?: "${projectDir}/assets/empty.bed", checkIfExists: true)
    exclusionMask = file(params.exclusion_bed ?: "${projectDir}/assets/empty.bed", checkIfExists: true)

    log.info """\
    ================================================================
    nf-FROG: Feedback-driven Reference Optimization through Genomic remapping
    samples                 : ${samples.size()}
    read units              : ${readUnits.size()}
    target chromosomes      : ${targetChromosomes.size()}
    maximum rounds          : ${maxRounds}
    minimum rounds          : ${minimumRounds}
    final heterozygotes     : ${finalHetMode}
    reference confidence    : ${gvcfMode}
    trim reads              : ${params.trim_reads}
    mark duplicates         : ${params.mark_duplicates}
    subset reference        : ${params.allow_reference_subset}
    ================================================================
    """.stripIndent()

    read_units_ch = Channel.fromList(readUnits)
        .map { row ->
            tuple(row.individual, row.unit, row.read1, row.read2, row.library,
                  row.karyotype, row.default_ploidy, row.pcr_free)
        }
    PREPROCESS_READ_UNIT(read_units_ch)

    if (params.kraken2_db) {
        kraken_input_ch = PREPROCESS_READ_UNIT.out
            .map { individual, unit, read1, read2, library, karyotype, default_ploidy, pcr_free, json, html ->
                tuple(individual, unit, read1, read2)
            }
        KRAKEN_SCREEN(kraken_input_ch, file(params.kraken2_db, checkIfExists: true))
    }

    SUBSET_REFERENCE(referenceFile, targetChromosomes)

    initial_refs_ch = Channel.fromList(samples)
        .combine(SUBSET_REFERENCE.out)
        .map { sample, target_reference ->
            tuple(target_reference.getName(), sample.individual, target_reference, target_reference)
        }
    targets_ch = Channel.fromList(callTargets)
        .map { target ->
            tuple(target.individual, target.chrom_order, target.chromosome,
                  target.shard_order, target.interval, target.ploidy)
        }
    sample_meta_ch = Channel.fromList(samples)
        .map { sample ->
            tuple(sample.individual, sample.karyotype, sample.default_ploidy, sample.pcr_free)
        }
    chromosome_meta_ch = Channel.fromList(
        targetChromosomes.withIndex().collect { targetChromosome, index ->
            tuple(String.format('%06d', index), targetChromosome)
        }
    )
    ploidy_mask_defs_ch = Channel.fromList(ploidyMaskDefs)
        .map { definition ->
            tuple(definition.individual, definition.chrom_order,
                  definition.chromosome, definition.intervals)
        }
    MAKE_PLOIDY_MASK(ploidy_mask_defs_ch)

    RUN_ITERATIONS(
        PREPROCESS_READ_UNIT.out,
        initial_refs_ch,
        targets_ch,
        sample_meta_ch,
        chromosome_meta_ch,
        MAKE_PLOIDY_MASK.out,
        1,
        maxRounds,
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.value(repeatMask),
        Channel.value(mappabilityMask),
        Channel.value(exclusionMask)
    )

    final_chromosome_input_ch = RUN_ITERATIONS.out.final_artifacts
        .map { individual, rnd, chrom_order, chromosome,
               final_vcf, final_tbi, final_mask, reason_mask,
               gvcf, gvcf_tbi, genotyped_vcf, genotyped_tbi,
               reference, fai, dict, filter_stats ->
            tuple(individual, rnd, chrom_order, chromosome,
                  final_vcf, final_tbi, final_mask, reference, fai)
        }
    BUILD_FINAL_CHROMOSOME(final_chromosome_input_ch)

    final_reference_ch = BUILD_FINAL_CHROMOSOME.out
        .groupTuple(by: 0)
        .map { individual, chrom_orders, chromosomes, fastas ->
            tuple(individual, fastas.flatten())
        }
        | BUILD_FINAL_REFERENCE

    final_vcfs_ch = RUN_ITERATIONS.out.final_artifacts
        .groupTuple(by: 0)
        .map { individual, rounds, chrom_orders, chromosomes,
               final_vcfs, final_tbis, final_masks, reason_masks,
               gvcfs, gvcf_tbis, genotyped_vcfs, genotyped_tbis,
               references, fais, dicts, filter_stats ->
            tuple(individual, final_vcfs.flatten())
        }
    GATHER_FINAL_VCF(final_vcfs_ch)

    final_gvcfs_ch = RUN_ITERATIONS.out.final_artifacts
        .groupTuple(by: 0)
        .map { individual, rounds, chrom_orders, chromosomes,
               final_vcfs, final_tbis, final_masks, reason_masks,
               gvcfs, gvcf_tbis, genotyped_vcfs, genotyped_tbis,
               references, fais, dicts, filter_stats ->
            tuple(individual, gvcfs.flatten())
        }
    GATHER_FINAL_GVCF(final_gvcfs_ch)

    final_masks_ch = RUN_ITERATIONS.out.final_artifacts
        .groupTuple(by: 0)
        .map { individual, rounds, chrom_orders, chromosomes,
               final_vcfs, final_tbis, final_masks, reason_masks,
               gvcfs, gvcf_tbis, genotyped_vcfs, genotyped_tbis,
               references, fais, dicts, filter_stats ->
            tuple(individual, final_masks.flatten(), reason_masks.flatten())
        }
    GATHER_FINAL_MASK(final_masks_ch)

    callable_input_ch = final_reference_ch
        .combine(GATHER_FINAL_MASK.out, by: 0)
        .map { individual, consensus, fai, noncallable_mask, reason_mask ->
            tuple(individual, consensus, fai, noncallable_mask)
        }
    MAKE_CALLABLE_BED(callable_input_ch)

    original_diff_input_ch = final_reference_ch
        .combine(SUBSET_REFERENCE.out)
        .map { individual, consensus, consensus_fai, original_reference ->
            tuple(individual, consensus, consensus_fai, original_reference)
        }
    ORIGINAL_REFERENCE_DIFF(original_diff_input_ch)

    final_qc_input_ch = final_reference_ch
        .combine(GATHER_FINAL_MASK.out, by: 0)
        .combine(GATHER_FINAL_VCF.out, by: 0)
        .map { individual, consensus, fai, mask, reason_mask, vcf, vcf_tbi ->
            tuple(individual, consensus, fai, mask, vcf)
        }
    FINAL_SAMPLE_QC(final_qc_input_ch)

    RUN_ITERATIONS.out.convergence_metrics.collect() | CONVERGENCE_SUMMARY
}
