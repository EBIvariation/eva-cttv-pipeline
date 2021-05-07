# Repeat expansion variants pipeline for ClinVar data

## 1. Introduction
This pipeline is a companion to the main [functional consequence & gene mapping pipeline](../vep_mapping_pipeline) and is intended to be used in conjunction with it. The main pipeline uses Ensembl's Variant Effect Predictor (VEP) to assign functional consequences to variants. The problem with VEP is that it does not support _repeat expansion variants._ They are a special class of duplication variants which represent uncertain number of repeats of a small set of nucleotides (see [Table 1 in this abstract](https://www.annualreviews.org/doi/10.1146/annurev.neuro.29.051605.113042) for examples). As a result of this limitation, the main pipeline cannot process these variants and ignores them. This sub-module was created to address this problem.

## 2. Running the pipeline
The pipeline has one input file and two output files.

### Input file
The input file is ClinVar's [compressed XML data dump](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz).

### Consequences table
The first output file, in TSV format, is the consequences table. It uses the same six-column format as the main pipeline. Example (tabs are replaced with spaces for readability):
```
RCV000986111    1    ENSG00000177570    SAMD12    short_tandem_repeat_expansion     0
RCV001003411    1    ENSG00000197386    HTT       trinucleotide_repeat_expansion    0
```

The columns, in order, are:
1. **RCV ID:** identifier of a ClinVar record containing the repeat variant. RCV IDs are used instead of the usual `CHROM:POS:REF:ALT` notation because it is not easily applied to repeat expansion variants.
2. Unused, retained for compatibility with the original pipeline.
3. **Ensembl gene ID.**
4. **Ensembl gene name.**
5. **Repeat expansion variant type,** can be either `short_tandem_repeat_expansion` or `trinucleotide_repeat_expansion`. See more on that below.
6. Unused, retained for compatibility with the original pipeline.

### Dataframe table
The second output file, also in TSV format, is the dump of the Pandas dataframe used for data processing. It contains a number of intermediate columns which were used to make the decision on the final data appearing in the consequences table. This table can be used for debugging or discussion purposes.

### Running the pipeline
See instructions [in the main README file](../README.md) to install the necessary packages first. Then run the pipeline:
```bash
python3 run_repeat_expansion_variants.py \
  --clinvar-summary-tsv variant_summary.txt.gz \
  --output-consequences output_consequences.tsv \
  --output-dataframe output_dataframe.tsv
```

## 3. Overview of data processing workflow
The information which the pipeline needs to determine for each variant is fairly simple:
* Whether the variant is a _trinucleotide_ expansion, or any other expansion, e. g. with the repeat unit length of 2.
* Which gene the variant affects.

This is made quite complicated by various systemic peculiarities in ClinVar data. In this section, the pipeline operation is explained step by step along with how those peculiarities are addressed.

### Step 1. Load and preprocess data
Repeat expansion variants are represented as **Microsatellite** events in ClinVar data. The first level of processing is to classify them according to their attributes:
* If the Microsatellite event has explicit coordinates and allelic sequences (CHROM/POS/REF/ALT), we can determine its type easily:
  + If len(REF) > len(ALT), this is a deletion (repeat contraction) event, which we're not interested in.
  + If len(REF) < len(ALT), but the difference is less than a certain threshold (see the current value in code as `REPEAT_EXPANSION_THRESHOLD`), it will be treated as a normal short insertion and not processed by this pipeline. This is done in this way to filter out a few thousands of events which constitute very short repeat expansions and are unlikely to be clinically relevant.
  + If the event passes the threshold, it is considered a repeat expansion variant with an explicit sequence. For it, we can easily calculate the length and see if it is a trinucleotide expansion or not.
* A small number of Microsatellite records do not have complete coordinates, and we have to guess their type by parsing their name, which is usually a HGVS-like expression. The rest of the pipeline, described below, deals mostly with those cases.

As input data from Clinvar, we get four useful columns:
* **Name:** variant identifier, which can be in three possible formats, explained below.
* **RCVaccession:** accession of a ClinVar record associating this variant with a phenotype.
* **GeneSymbol:** symbol (name) of a gene which the variant affects.
* **HGNC_ID:** HGNC ID of a gene which the variant affects

**RCVaccession** and **GeneSymbol** columns may contain multiple values, delimited by semicolons. This means there are multiple records associating the variant with a phenotype, or the variant affecting multiple genes. For example, variant `NM_001007026.1(ATN1):c.1462CAG[(49_55)] (p.Gln488[(49-55)])` is associated with accessions RCV000032095, RCV000032096, RCV000032097 and gene symbols ATN1 and LOC109461484.

In case there is a single gene in **GeneSymbol** column, then **HGNC_ID** will contain the ID of that gene. If there are multiple genes, the **HGNC_ID** column will be empty.

As the first preprocessing step, we split each record by **RCVaccession** and **GeneSymbol** columns. This means that, in the example above, a single record becomes six, one per every combination of ClinVar record and gene.

Next, we remove duplicates in the resulting table. They can occur because some (Name, RCVaccession, GeneSymbol) tuples appeared twice in the original table, with coordinates for GRCh37 and GRCh38. Since the genomic coordinates are not used by this pipeline, we simply remove the duplicates.

### Step 2. Parse variant identifiers
The **Name** column contains a variant description which can be in three possible formats:
1. HGVS-like variant description in coding or genomic coordinates, e. g. `NM_000044.4(AR):c.172_174CAG(10_36) (p.Gln69_Gln80del)`.
2. HGVS-like variant description in protein coordinates, e. g. `NP_003915.2:p.Ala260(5_9)`.
3. Human-readable, however standardised, description, e. g. `ATN1, (CAG)n REPEAT EXPANSION`.

There are also two properties which we would like to extract:
1. Length of the repeat unit. It can be determined using two approaches:
  + Directly from the sequence. In examples 1 and 3 it would be len(CAG) = 3.
  + Indirectly from the coordinate span of the variant. In example 1, it would be 174 – 172 + 1 = 3.
2. Transcript in which the variant resides. In example 1, it would be NM_000044.

Here is the table describing which properties can be extracted from which identifier types:

| Identifier type           | Repeat sequence | Coordinate span | Transcript ID |
| ------------------------- | --------------- | --------------- | ------------- |
| HGVS-like, coding/genomic | sometimes*      | sometimes*      | yes           |
| HGVS-like, protein        | –               | –               | –             |
| Free text                 | yes             | –               | –             |

_\* At least one of (repeat sequence, coordinate span) is always present in this case._

At this step, we extract all properties which we can from the variant identifiers and save them to separate dataframe columns.

### Step 3. Annotate Ensembl gene information
There are three columns which can be used to deduce gene information: HGNC_ID & GeneSymbol (from ClinVar data), and TranscriptID (parsed from the variant identifier). Because neither of them is present for **all** records, we successively try to annotate Ensembl gene ID based on each of those columns in that order. For example, if Ensembl gene ID for a particular record was determined based on HGNC_ID, then GeneSymbol and TranscriptID will not be tried for that column. Ensembl gene information is being annotated using their [BioMart service](https://www.ensembl.org/biomart/).

In some cases, a mapping from HGNC ID produces multiple Ensembl hits. For example, HGNC:10560 is being resolved to ENSG00000285258 and ENSG00000163635. In this case, a record is further split into two, each mapping into their corresponding Ensembl gene IDs.

### Step 4. Determine repeat type and record completeness
There are two repeat types detected by this pipeline:
* `trinucleotide_repeat_expansion`, corresponding to SO term http://sequenceontology.org/browser/current_release/term/SO:0002165. This is chosen when the repeat length appears to be either 3, **or a multiple or 3.** For example:
  + `NM_001081563.1:c.*224_*226[(50-?)]` — repeat unit spans 3 nucleotides based on coordinates
  + `NM_001081560.2(DMPK):c.*224_*226CTG(51_?)` — spans 3 nucleotides based on coordinates, and the sequence length is also 3
  + `NM_000100.3(CSTB):c.-210CCCCGCCCCGCG(2_3)` — sequence length is 12, a multiple of 3
* `short_tandem_repeat_expansion`, corresponding to SO term http://sequenceontology.org/browser/current_release/term/SO:0002162. This is chosen in the cases when the repeat length can be _determined_ by at least one method, but its unit length is not a multiple of 3.

Similarly to the previous step, repeat type can be determined using several fields. The determination algorithm is as follows:

* If the variant identifier is type 3, “Protein HGVS”: assume repeat is a `trinucleotide_repeat_expansion`, since it affects entire amino acids. (Each amino acid is encoded by a triplet of nucleotides. So if we see an insertion of N amino acids, it means an insertion of 3×N nucleotides took place, hence this will be a trinucleotide repeat expansion.)
* Otherwise:
  + Determine repeat unit length
    - If available, use length determined directly from sequence as a priority
    - Otherwise, use coordinate span
  + If repeat unit length is a multiple of 3, assign `trinucleotide_repeat_expansion`. Otherwise, assign `short_tandem_repeat_expansion`.

In case there is no evidence for either repeat type, the field will be left empty.

The record defined as being complete if all three fields are present: variant type; Ensembl gene ID; Ensembl gene name.

### Step 5. Generating the output files
Dataframe is output as is to faciliate debugging and discussion. It includes both complete and incomplete records. Example of such a dataframe, generated on 2020-04-08: https://docs.google.com/spreadsheets/d/19rDkxD8rEbVKhIpg-1EyHuTSeFP4IrSmcm7PkO1bLIg/edit?usp=sharing.

The consequences table is obtained by removing the unnecessary columns and collapsing the dataframe. The dataframe is grouped by all triples of (RCVaccession, EnsemblGeneID, EnsemblGeneName), and there is a check to ensure that for each triple there is only one type of repeat predicted.

The consequences table only includes complete records and is saved in a six-column format suitable for use in the main Open Targets evidence string generation pipeline.

## Note on incomplete records
As of 2020-04-08, the pipeline was able to process 111 out of 115 records. The causes for missing the 4 records are:
* ClinVar data issues
  + RCV000001021 has non-standard name “fragile site, folic acid type, rare, fra(12)(q13.1)” which cannot be parsed.
  + RCV000192065 has an empty name.
* Ensembl missing a gene: RCV000000215 and RCV000006519 are mapped to two genes each, `ATXN8` and `ATXN8OS`, generating four RCV/gene record pairs. Mappings to `ATXN8OS` are processed normally; however, the `ATXN8` gene is missing from Ensembl 99 release.
