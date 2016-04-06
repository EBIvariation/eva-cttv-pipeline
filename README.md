## README ##

[![Build Status](https://travis-ci.org/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.org/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)


Minimum Python version needed: 3.4


Building and (optional) Setting up virtual environment
-------

1. "git clone --recursive git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. [OPTIONAL] "virtualenv -p python3.4 venv"
4. [OPTIONAL] "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
5. pip install -r requirements.txt
6. and then:
	7. to install: "python3 setup.py install"
	8. to install to develop: "python3 setup.py develop"
	9. to build a source distribution: "python3 setup.py sdist"


Usage
-------

python3 bin/main.py --out \<OUTPUT_FILE\> -e \<EFO_MAPPING_FILE\> -g \<SNP_2_GENE_MAPPING_FILE\> -v \<VARIANT_SUMMARY_FILE\> [--clinSig \<CLINICAL_SIGNIFICANCE_LIST\>] [--ignore \<TERM_URL_IGNORE_FILE\>]

OUTPUT_FILE: path to output directory

EFO_MAPPING_FILE: file with mappings of Clinvar trait names to URLS (e.g. resources/ClinVar_Traits_EFO_090915.xls)

SNP_2_GENE_MAPPING_FILE: file from CTTV with mappings from rs IDs to Ensembl gene IDs, and to functional consequences (e.g. resources/cttv012_snp2gene_20160222.tsv)

VARIANT_SUMMARY_FILE: variant_summary file from Clinvar (e.g. resources/variant_summary_2015-05.txt)

CLINICAL_SIGNIFICANCE_LIST: comma separated (no spaces) list of clinical significances allowed to generate evidence strings (defaults to "pathogenic,likely pathogenic")

TERM_URL_IGNORE_FILE: path to file containing list of invalid EFO URLs which will be used to filter out evidence strings with matching URLs; file should contain list of URLs, one on each line
