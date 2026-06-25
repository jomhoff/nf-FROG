#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# nf-FROG uses recursive workflow composition and requires the legacy
# parser. Nextflow 26 defaults to a parser that does not compile this pattern.
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

exec nextflow run "${script_dir}/main.nf" \
  -c "${script_dir}/nextflow.config" \
  "$@"
