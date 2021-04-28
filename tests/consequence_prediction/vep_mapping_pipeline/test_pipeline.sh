#!/bin/bash

# TODO fix this

set -Eeuo pipefail

# VEP / repeat expansion pipeline tests.
# For the actual test, we're running a set of 2,000 ClinVar variants through VEP and comparing the result with the
# expected one (diff will exit with return code 0 if the files are identical, and with 1 otherwise). Of course, this
# means that when VEP updates, the test will break; however, this is exactly the intention, as in this case we will be
# able to compare the results and see if they make sense.
echo 'Testing VEP mapping pipeline'
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $SCRIPT_DIR
#SOURCE_DIR="$(dirname $(dirname $SCRIPT_DIR))/nextflow"
bash run_consequence_mapping.sh vep_mapping_pipeline/test/input.xml.gz output_mappings.tsv
diff vep_mapping_pipeline/test/output_mappings.tsv output_mappings.tsv