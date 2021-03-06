# Manual curation, part III, technical: export curation results and submit feedback to ZOOMA

Before running, set up the environment:
* [Common environment](../environment.md)
* [Protocol-specific environment](README.md#setting-up-environment)

## Extract curation results from the spreadsheet

### Finished mappings
Once the manual curation is completed, apply a spreadsheet filter so that only traits with **Status = DONE** are visible. Copy data for all non-empty rows from three columns: “ClinVar label”; “URI of selected mapping”; “Label of selected mapping”, in that order. **Do not include header lines.** Save the data to a file `${CURATION_RELEASE_ROOT}/finished_mappings_curation.tsv`.

### Terms requiring import into EFO
After the manual curation has been completed, traits remaining unmapped or poorly mapped should be submitted to EFO if a suitable parent term is available. Open the curation spreadsheet and use filters to display only terms with **Status = IMPORT.** Copy _just the ontology_ URLs into the file `${CURATION_RELEASE_ROOT}/terms_for_efo_import.txt`, one URL per line, for example:
```
http://purl.obolibrary.org/obo/HP_0002647
http://purl.obolibrary.org/obo/MONDO_0000727
http://www.orpha.net/ORDO/Orphanet_199306
```

## Run the automated protocol

```bash
# Define some variables for readability
export EXISTING_MAPPINGS=${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv
export NEW_MAPPINGS=${CURATION_RELEASE_ROOT}/trait_names_to_ontology_mappings.tsv

# Concatenate finished automated and manual mappings into a single file
cat \
  ${CURATION_RELEASE_ROOT}/automated_trait_mappings.tsv \
  ${CURATION_RELEASE_ROOT}/finished_mappings_curation.tsv \
| sort -u > ${NEW_MAPPINGS}

# Add all mappings from the database which are *not* present in the results of the current curation iteration (automated
# + manually curated). This is done in order to never lose mappings, even if they are not present in ClinVar during the
# latest curation iteration.
# The first file operand is the list of mappings in the current database; and the second is the list of trait names
# which are only present in the existing database and not in the new mappings.
export LC_ALL=C
join -j 1 -t $'\t' \
  <(sort -t $'\t' -k 1,1 ${EXISTING_MAPPINGS}) \
  <(comm -23 <(cut -d $'\t' -f 1 ${EXISTING_MAPPINGS} | sort -u) <(cut -d $'\t' -f 1 ${NEW_MAPPINGS} | sort -u)) \
>> ${NEW_MAPPINGS}

# Run the helper script to prepare the table for EFO import
python3 ${CODE_ROOT}/bin/trait_mapping/create_efo_table.py \
  -i ${CURATION_RELEASE_ROOT}/terms_for_efo_import.txt \
  -o ${CURATION_RELEASE_ROOT}/efo_import_table.tsv

# Generate ZOOMA feedback
echo -e 'STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE' \
  > ${CURATION_RELEASE_ROOT}/eva_clinvar.txt
tail -n+2 ${NEW_MAPPINGS} \
  | cut -f-2 \
  | sort -t$'\t' -k1,1 \
  | awk -F$'\t' -vDATE="$(date +'%y/%m/%d %H:%M')" '{print "\t\tdisease\t" $1 "\t" $2 "\teva\t" DATE}' \
>> ${CURATION_RELEASE_ROOT}/eva_clinvar.txt

# Update the symbolic links pointing to the location of the most recent curation result and ZOOMA feedback dataset.
# This will be used by the main evidence string generation protocol.
ln -s -f ${NEW_MAPPINGS} ${EXISTING_MAPPINGS}
ln -s -f ${CURATION_RELEASE_ROOT}/eva_clinvar.txt ${BATCH_ROOT_BASE}/manual_curation/eva_clinvar.txt
```

## Check that the resulting file contains no duplicates
The resulting list of text-to-ontology mappings should not contain any complete duplicates. Check that this is the case by using the following command. If everything is correct, it should not output anything:

```bash
sort ${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv | uniq -c | awk '$1 > 1'
```

If there are duplicates, resolve this by editing the `${BATCH_ROOT_BASE}/manual_curation/latest_mappings.tsv` file directly.

## Copy the table for EFO import
The file `${CURATION_RELEASE_ROOT}/efo_import_table.tsv` will contain a partially ready table for EFO import. Copy its contents into the “Add EFO disease” sheet in the curation spreadsheet.

## Submit feedback to ZOOMA
See more details on ZOOMA feedback in the [evidence string generation protocol](../generate-evidence-strings.md#submit-feedback-to-zooma). At this stage, only the **eva_clinvar** dataset is being submitted; clinvar_xrefs is submitted during evidence string generation.

```bash
# EXECUTE UNDER FTP ADMINISTRATIVE USER
# DON'T FORGET TO SET THE TWO VARIABLES BELOW AGAIN
export BATCH_ROOT_BASE=...
export FTP_PATH_BASE=...

# Create the folder, copy the file to FTP, and update the “latest” folder
FTP_PATH=${FTP_PATH_BASE}/`date +%Y/%m/%d`
mkdir -p ${FTP_PATH}
cp ${BATCH_ROOT_BASE}/manual_curation/eva_clinvar.txt ${FTP_PATH}
cp ${FTP_PATH}/eva_clinvar.txt ${FTP_PATH_BASE}/latest/eva_clinvar.txt
```

After uploading both files, confirm that the changes have propagated to the FTP:
```bash
md5sum ${BATCH_ROOT_BASE}/manual_curation/eva_clinvar.txt
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/`date +%Y/%m/%d`/eva_clinvar.txt | md5sum
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt | md5sum
```
