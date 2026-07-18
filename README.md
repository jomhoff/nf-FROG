# nf-FROG

**nf-FROG** — **F**eedback-driven **R**eference **O**ptimization through
**G**enomic remapping — builds callable, individual-specific genome consensus
sequences from paired-end short reads through iterative reference refinement.

It maps each sample to the current personalized reference, calls variants and
reference-confidence blocks with FreeBayes, applies ploidy-aware genotype and variant
filters, updates the mapping reference with confident homozygous SNPs, and
repeats until the reference converges. A separate final step masks every site
that cannot be defended as callable and represents passing heterozygous SNPs as
IUPAC codes, a major allele, or `N`.

## Why nf-FROG?

Ordinary variant-only consensus generation leaves unchanged reference bases at
sites where no trustworthy reference call was made. That can make missing or
uncallable sequence look like confidently observed reference sequence.

`nf-FROG` uses the gVCF reference-confidence model to distinguish:

- confident reference calls;
- confident alternate calls;
- heterozygous calls; and
- sites that must be masked because evidence is missing or unreliable.

The iterative mapping reference and final masked consensus are intentionally
different products. Intermediate references are never masked, so a low-quality
site in one round remains recoverable in a later round.

## Main features

- Multiple samples, lanes, and libraries
- Per-read-unit read groups
- Optional read trimming with `fastp`
- Duplicate marking with `samtools markdup`
- Optional Kraken2 contamination screening
- Ploidy-aware calling by sample, chromosome, and interval
- Explicit support for ploidy-zero intervals
- Parallel per-interval FreeBayes calls with gVCF reference blocks
- MNP, complex-call, and multiallelic decomposition with `bcftools norm`
- Robust chromosome-specific depth thresholds recalculated each round
- Site, genotype, allele-depth, and allele-balance filtering
- Reason-coded noncallable masks
- Optional repeat, low-mappability, and exclusion masks
- Sample-specific convergence and oscillation detection
- Original-reference-coordinate SNP VCF reconstruction
- Mapping, depth, insert-size, duplication, Ti/Tv, allele-balance, and
  callability QC

## Biological model

### Iterative references

Only passing homozygous alternate SNPs change the iterative reference.
Heterozygous calls do not alter it, and uncertain sites retain the current
reference base. This keeps mapping references in `A/C/G/T/N` space and avoids
random allele-depth fluctuations causing heterozygous sites to flip between
rounds.

### Final consensus

The final consensus is made from the reference used in the terminal calling
round and the terminal filtered VCF. Passing VCF alleles are applied before the
final mask, preserving VCF REF/FASTA agreement while still turning noncallable
sites into `N`. Its noncallable mask is the union of:

- gVCF-absent spans;
- no-call or missing-depth records;
- exact low- and high-depth positions from `samtools depth`;
- rejected variant sites;
- ploidy-zero intervals;
- pre-existing non-ACGT reference sequence;
- optional repeat masks;
- optional low-mappability masks; and
- optional user exclusion masks.

A rejected alternate call is therefore masked. It never silently becomes a
confident reference base.

### Coordinates and indels

The FASTA products remain SNP-only so their coordinates match the supplied
reference. Indels and complex variants remain available in the raw FreeBayes VCF
and gVCF, but they are not applied to the consensus FASTA.

### Heterozygous SNPs

`--final_het_mode` controls passing heterozygous SNPs:

| Mode | Result |
| --- | --- |
| `iupac` | Genotype-based IUPAC ambiguity code; default |
| `major` | Allele with greatest `AD`; REF wins ties |
| `mask` | `N` |

This is a consensus representation, not a phased diploid assembly. It cannot
represent phase, heterozygous indels, or structural variants.

## Requirements

- Nextflow 24.10–25.x
- Java 17 or newer
- Conda/Mamba for the supplied profile, or equivalent tools installed manually

The provided environment includes:

- FreeBayes
- bwa-mem2
- samtools
- bcftools/htslib
- bedtools
- fastp
- Kraken2
- Python

The workflow currently uses Nextflow's recursion-capable `v1` syntax parser.
Use the supplied launcher, which sets this automatically.

## Installation

Clone or copy the `nf-FROG` directory, then make the launcher executable:

```bash
chmod +x run.sh
```

The Conda environment is created automatically when running with
`-profile conda`.

## Inputs

### Reference FASTA

Provide the starting reference with `--ref_file`.

Sequence identifiers must contain only letters, numbers, periods, underscores,
and hyphens. Reference order is preserved in the output.

### Chromosome list

`--chromos_file` is a plain-text file containing one reference sequence name per
line:

```text
chr1
chr2
chrZ
chrW
```

By default, every reference contig must be listed. To intentionally analyze and
emit only a subset, use `--allow_reference_subset true`.

### Sample sheet

A tab-separated sample sheet is recommended:

```text
individual	read1	read2	read_group	library	karyotype	default_ploidy	pcr_free
bird01	/data/bird01_L1_R1.fastq.gz	/data/bird01_L1_R2.fastq.gz	bird01_L1	lib1	ZZ	2	true
bird01	/data/bird01_L2_R1.fastq.gz	/data/bird01_L2_R2.fastq.gz	bird01_L2	lib1	ZZ	2	true
bird02	/data/bird02_R1.fastq.gz	/data/bird02_R2.fastq.gz	bird02_L1	lib2	ZW	2	true
```

Required columns:

| Column | Meaning |
| --- | --- |
| `individual` | Biological sample identifier |
| `read1` | R1 FASTQ/FASTQ.GZ |
| `read2` | R2 FASTQ/FASTQ.GZ |

Optional columns:

| Column | Default | Meaning |
| --- | --- | --- |
| `read_group` | Generated | Unique sequencing unit identifier |
| `library` | Individual ID | Library identifier used for duplicate marking |
| `karyotype` | `default` | Key used by the ploidy table |
| `default_ploidy` | `2` | Ploidy where no explicit interval rule exists |
| `pcr_free` | `false` | Retained for sample-sheet compatibility; duplicate marking is controlled globally |

`individual`, `read_group`, and `library` may contain only letters, numbers,
periods, underscores, and hyphens. Read-group values must be globally unique.

The legacy `--indivs_file` plus `--readsdir` interface is retained, but cannot
describe multiple units, libraries, karyotypes, or PCR-free status.

### Ploidy table

Use `--ploidy_file` to define sex chromosomes, pseudoautosomal regions,
haploid regions, or absent chromosomes:

```text
karyotype	chromosome	start	end	ploidy
ZZ	chrZ	1	*	2
ZZ	chrW	1	*	0
ZW	chrZ	1	500000	2
ZW	chrZ	500001	*	1
ZW	chrW	1	*	1
```

Coordinates are one-based and inclusive. `*` means the chromosome end.

For a karyotype/chromosome with explicit rows, the rows must form a complete,
gap-free, non-overlapping partition. Ploidy-zero intervals are fully masked and
cannot alter the iterative reference.

### Optional masks

BED3 or BED4 files may be supplied with:

- `--repeat_mask_bed`
- `--low_mappability_bed`
- `--exclusion_bed`

Coordinates must match the starting reference.

## Quick start

Local execution:

```bash
./run.sh \
  -profile conda,local \
  --samplesheet samples.tsv \
  --ploidy_file ploidy.tsv \
  --chromos_file chromosomes.txt \
  --ref_file reference.fa \
  --outputdir results_nf_FROG
```

SLURM execution:

```bash
./run.sh \
  -profile conda,slurm \
  --samplesheet samples.tsv \
  --ploidy_file ploidy.tsv \
  --chromos_file chromosomes.txt \
  --ref_file reference.fa \
  --outputdir results_nf_FROG
```

Resume an interrupted run with `-resume`.

## Important parameters

### Iteration and convergence

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--max_rounds` | 6 | Maximum mapping/calling rounds |
| `--minimum_rounds` | 2 | Earliest round eligible for convergence |
| `--convergence_max_changes` | 0 | Maximum changed bases at convergence |

A sample converges when the change threshold is met and no two-cycle
oscillation is detected. Samples retire independently.

### Mapping resources

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--samtools_sort_memory` | `1G` | Memory per samtools sort thread during read-unit BAM generation |

The default SLURM mapping profile requests 8 CPUs and 96 GB RAM. Large
references can exceed the earlier 16 CPU / 48 GB profile because `samtools sort`
allocates memory per thread while `bwa-mem2` is also resident. If mapping jobs
exit with status 137 or show `samtools sort ... Killed`, increase the mapping
process memory in a site config and/or lower `--samtools_sort_memory`.

### Calling and filtering

| Parameter | Default |
| --- | ---: |
| `--freebayes_region_size` | 10,000,000 |
| `--freebayes_gvcf_chunk` | 50 |
| `--freebayes_min_alternate_count` | 2 |
| `--freebayes_min_alternate_fraction` | 0.05 |
| `--freebayes_use_best_n_alleles` | 4 |
| `--min_qual` | 30 |
| `--min_mqm` | 20 |
| `--min_alt_mean_quality` | 20 |
| `--balanced_alt_count` | 4 |
| `--min_alt_strand_observations` | 1 |
| `--min_alt_placement_observations` | 1 |
| `--min_het_allele_fraction` | 0.20 |
| `--max_het_major_fraction` | 0.80 |
| `--min_hom_alt_fraction` | 0.90 |
| `--min_allele_depth` | 3 |

These are starting points, not organism-independent truth. Calibrate them for
coverage, divergence, library construction, assembly repetitiveness, expected
heterozygosity, and available validation data.

`--freebayes_region_size` controls SLURM parallelism. Each chromosome or
scaffold is divided into non-overlapping regions of at most this size while
respecting ploidy boundaries. FreeBayes itself uses one CPU per region, allowing
Nextflow and SLURM to distribute regions independently.

### Depth and callability

| Parameter | Default |
| --- | ---: |
| `--absolute_min_depth` | 5 |
| `--min_depth_multiplier` | 0.5 |
| `--max_depth_multiplier` | 2.0 |
| `--depth_mad_multiplier` | 6 |
| `--depth_upper_quantile` | 0.99 |
| `--mapq_min` | 20 |
| `--baseq_min` | 20 |

Thresholds are chromosome-specific and recalculated each round from positive
depth observations. Sparse chromosomes use a genome-wide fallback.

FreeBayes gVCF blocks establish whether the caller emitted evidence across each
reference span. The default 50-base block size keeps files manageable. Exact
low/high-depth masking is calculated independently at every base with
`samtools depth`, so block compression cannot hide local coverage failures.
Final mask merging uses a disk-backed sort/sweep so chromosomes with many
small rejected or low-depth intervals do not require keeping all intervals in
Python memory.

### Preprocessing

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--trim_reads` | `false` | Enable fastp trimming; QC always runs |
| `--mark_duplicates` | `true` | Mark PCR/optical duplicates |
| `--kraken2_db` | unset | Run optional contamination screening |

## Output layout

```text
results_nf_FROG/
├── 00.reference/
├── 00.reference_indexes/
├── 00.read_qc/
├── 01.bams/
├── 02.gvcfs/
├── 02.variants/
├── 03.consensus/
├── 04.consensus_refs/
├── 05.round_consensus/
├── 06.coverage/
├── 07.qc/
└── 08.masks/
```

Principal per-sample outputs under `03.consensus/<sample>/` include:

| Output | Meaning |
| --- | --- |
| `<sample>.consensus.fa` | Final masked consensus |
| `<sample>.final.g.vcf.gz` | Terminal gVCF |
| `<sample>.final.current_reference.vcf.gz` | Terminal filtered VCF |
| `<sample>.original_reference.vcf.gz` | SNP differences reconstructed against the starting reference |
| `<sample>.callable.bed` | Callable intervals |
| `<sample>.noncallable_mask.bed` | Union mask |
| `<sample>.noncallable_mask.reasons.bed` | Mask with failure reasons |

Per-round outputs include BAMs, FreeBayes gVCFs, raw decomposed VCFs, iterative
references, depth thresholds, filter metrics, reason-specific masks, SAMtools
alignment and insert-size metrics, variant QC, and convergence statistics.

## QC interpretation

Inspect at least:

- callable fraction;
- mapping rate and mismatch rate across rounds;
- duplicate fraction;
- insert-size distributions;
- chromosome-specific depth thresholds;
- changed and oscillating bases per iteration;
- heterozygous reference-allele fraction;
- Ti/Tv;
- filter-reason counts; and
- mask-reason totals.

Large round-to-round changes, persistent oscillations, strong heterozygous
reference imbalance, unusual Ti/Tv, or extensive high-depth masking can signal
reference divergence, paralogy, contamination, ploidy errors, or unsuitable
hard filters.

## Known limitations

- Designed for germline paired-end short-read data
- Not a mitochondrial heteroplasmy workflow
- Does not produce phased haplotypes
- Does not apply indels or structural variants to FASTA
- Does not correct assembly structure or collapsed paralogs
- Hard-filter defaults require study-specific validation
- Recursive workflow composition is a Nextflow preview feature

For mitochondrial analysis, use a workflow designed for circular genomes,
NUMTs, and low-frequency heteroplasmy.

## Testing

Run helper tests with:

```bash
python3 -m unittest tests/test_helpers.py
```

Generate the small deterministic end-to-end test genome with:

```bash
python3 tests/generate_synthetic_dataset.py \
  --output tests/synthetic_run/input
```

After running nf-FROG on that dataset, validate its biological outputs with:

```bash
python3 tests/validate_synthetic_run.py \
  --input tests/synthetic_run/input \
  --output tests/synthetic_run/output
```

Validate the workflow graph without executing tasks:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
  -c nextflow.config \
  -profile local \
  -preview \
  --samplesheet samples.tsv \
  --chromos_file chromosomes.txt \
  --ref_file reference.fa
```

## License and citation

If nf-FROG is used in published work, cite the versions of Nextflow, FreeBayes,
bwa-mem2, SAMtools,
BCFtools, BEDTools, and fastp used for the analysis.
