#!/bin/bash
set -euxo pipefail

# A wrapper script for running the consequence mapping pipeline concurrently, based on ClinVar XML input.
# See /README.md for detailed information.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# If you're running on VCF, substitute the first line with: bcftools query "$1" -f '%CHROM:%POS:%REF:%ALT\n'
python3 "${DIR}/extract_variants_for_vep.py" --clinvar-xml "$1" \
  | sort -u \
  | parallel \
    --halt now,fail=1    `# If any job fails, kill the remaining ones immediately and report failure`    \
    --pipe               `# Input is read from STDIN and split by chunks`                                \
    -j 20                `# Number of concurrent workers`                                                \
    -N 200               `# Number of records (lines) per worker`                                        \
    --tmpdir .           `# Store temporary files in the current directory to avoid /tmp overflow`       \
    python3 "${DIR}/vep_mapping_pipeline/consequence_mapping.py" \
  | sort -u > "$2"
