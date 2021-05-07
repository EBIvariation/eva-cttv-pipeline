# Gene & functional consequence mapping pipeline

## Outline

The pipeline located in the [vep_mapping_pipeline](/vep_mapping_pipeline) directory performs gene & functional consequence mapping. That is, it takes variants as input and maps each of them to a gene (Ensembl gene ID & gene name) and to the most severe functional consequence which a variant causes in that gene, according to Variant Effect Predictor.

## Running the pipeline, input & output formats

### Input file
The pipeline uses succinct, VCF-compatible variant identifiers, in the format of `CHROM:POS:REF:ALT`. Input file must contain variants in exactly this format, one entry per line. Example:
```
10:110821920:G:A
10:110821933:C:T
10:110821934:G:A
```

### Output file
Output is a TSV file consisting of six columns:
1. Variant identifier: the same one as in the input files.
2. The second column is not used and is always set to 1. It is retained for compatibility purposes (see more about that in the release notes).
3. Ensembl gene ID.
4. Ensembl gene name.
5. Most severe functional consequence of that variant for that gene.
6. Distance from a variant to a gene. The contents of this field will vary depending on the pipeline options:
  * By default, distance to the gene is not reported, and the field will always be zero.
  * With `--report-distance` flag, this will be nonzero only for upstream and downstream gene variants, and zero for all others.

In general, each variant can have more than one record in the output file, because it may affect multiple genes.

Example (note that tabs have been replaced with spaces here for readability):
```
10:110821920:G:A    1    ENSG00000203867    RBM20    missense_variant         0
10:110821933:C:T    1    ENSG00000203867    RBM20    missense_variant         0
10:110821934:G:A    1    ENSG00000203867    RBM20    splice_region_variant    0
```

### Running the script directly (on a small number of entries, for testing/debugging only)
The pipeline core module is [consequence_mapping.py](/vep_mapping_pipeline/consequence_mapping.py). It reads data from STDIN and writes to STDOUT in the formats described above. It can be run as follows:
```bash
python3 consequence_mapping.py <input_variants.txt >output_mappings.tsv
```

This should only be done for testing purposes and only for a small number of variants, because when querying VEP API, the script submits of all of the variants it receives in a single query (this is more efficient than submitting them one by one).

### Running the pipeline using a wrapper script
In production environment the pipeline should be run using a wrapper script which would take care of preprocessing and parallelisation. There is a simple wrapper script available, [run_consequence_mapping.sh](/run_consequence_mapping.sh). It can be run as follows:
```bash
bash run_consequence_mapping.sh input_variants.xml.gz output_mappings.tsv
```

It performs the following steps:
1. Take ClinVar's [compressed XML data dump](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz) as input.
1. Extract the necessary fields using `bcftools query`.
1. Sort and remove duplicates.
1. Using `parallel`, split the input file into chunks at most **200** records long and process them using **20** workers in parallel. These values have been selected to process the variants in reasonable time, but not to overload VEP API.
1. Sort the result and remove any duplicate values if present, and save it to a specified output file.

## Mapping process
For each variant:
1. Query VEP with default distance to the gene (5,000 bases either way). Consider only consequences which affect **either a protein coding transcript or a miRNA**. Find the most severe consequence for either of those transcripts, according to the list of priorities described [on Ensembl website](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html). Output **all non-redundant** consequences of this (most severe) type.
2. *Optionally:* if we found nothing during the previous step, repeat VEP search for this variant, now with a distance up to 500,000 bases either way. If there are any consequences which affect a **protein coding** transcript, choose the most severe consequence type (usually this will be either an upstream or a downstream gene variant) and output **a single consequence with the smallest distance**.

Distant querying (Step 2), while disabled by default, can be enabled by appending `--enable-distant-querying` flag to the core pipeline.

## Note on porting & changes from the original pipeline
This pipeline originated as a Python port of the original Open Targets ["SNP to gene"](https://github.com/opentargets/snp_to_gene) pipeline. An effort has been taken to retain backwards compatibility where possible; however, many important changes have been introduced. Please see the release notes for description of those changes.

## Note on external use
The pipeline is implemented to be as standalone as possible, in hope that it could in future become useful to other people (including, but not limited to, other Open Targets submitters). Feel free to contact the maintainers of this repository if you need any assistance in adopting this pipeline for your use case.
