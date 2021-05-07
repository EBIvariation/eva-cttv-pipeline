#!/bin/bash

set -Eeuo pipefail

# VEP pipeline tests.
# For the actual test, we're running a set of 2,000 ClinVar variants through VEP and comparing the result with the
# expected one (diff will exit with return code 0 if the files are identical, and with 1 otherwise). Of course, this
# means that when VEP updates, the test will break; however, this is exactly the intention, as in this case we will be
# able to compare the results and see if they make sense.
echo 'Testing VEP mapping pipeline'
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_DIR="$(dirname $(dirname $(dirname "${SCRIPT_DIR}")))"

bash "${ROOT_DIR}/consequence_prediction/run_consequence_mapping.sh" "${SCRIPT_DIR}/resources/input.xml.gz" output_mappings.tsv
diff "${SCRIPT_DIR}/resources/output_mappings.tsv" output_mappings.tsv

rm output_mappings.tsv