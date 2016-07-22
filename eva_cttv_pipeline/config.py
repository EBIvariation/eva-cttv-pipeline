LOCAL_SCHEMA = "resources/schema_local"


SPARQLEP = ""
SPARQL_USER = ""
SPARQL_PASSWORD = ""


#########################################
# clinvar_to_evidence_strings.py settings
#########################################
# cellbase
BATCH_SIZE = 200
HOST = 'wwwdev.ebi.ac.uk'

# output settings
EVIDENCE_STRINGS_FILE_NAME = 'evidence_strings.json'
EVIDENCE_RECORDS_FILE_NAME = 'evidence_records.tsv'
UNMAPPED_TRAITS_FILE_NAME = 'unmappedTraits.tsv'
UNAVAILABLE_EFO_FILE_NAME = 'unavailableefo.tsv'
NSV_LIST_FILE = 'nsvlist.txt'

######


##############################
# evidence_strings.py settings
##############################
GEN_EV_STRING_JSON = "resources/CTTVGeneticsEvidenceString.json"
SOM_EV_STRING_JSON = "resources/CTTVSomaticEvidenceString.json"
GEN_SCHEMA_FILE = LOCAL_SCHEMA + "/src/genetics.json"
SOM_SCHEMA_FILE = LOCAL_SCHEMA + "/src/literature_curated.json"

#############
