# Consequence prediction pipelines for ClinVar/Open Targets

This repository contains two modules and their corresponding wrapper scripts:
* [vep_mapping_pipeline](vep_mapping_pipeline): maps variants (CHROM:POS:REF:ALT) to their most severe functional
  consequence according to Ensembl VEP, as well as their Ensembl gene ID and name.
* [repeat_expansion_variants](repeat_expansion_variants): parses ClinVar `variant_summary` file and extracts
  information about repeat expansion variants.

Please see the corresponding module README file for more information.

## Installing requirements

The commands below has been tested for Ubuntu 18.04 and newer. You might have to adjust commands and package names if
you're using a different distribution. Note in particular that some older Debian and Ubuntu distributions include
ancient htslib/samtools/bcftools versions, which will not work.

```bash
sudo apt -y install samtools bcftools parallel libbz2-dev liblzma-dev
sudo python3 -m pip -q install -r requirements.txt
```