"""
Ensembl BioMart is an API which can be used to query Ensembl databases. It can be used to map external identifiers
(HGNC gene ID, gene symbol, RefSeq transcript ID) to Ensembl gene ID and name. The BioMart API is a bit quirky in that
it uses XML to specify the request.
"""

from io import StringIO
import logging

import pandas as pd
import requests
from retry import retry

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# The idea behind this template is that we query BioMart with a list of identifiers (`identifier_list`) from the
# `key_column`. For example, it can be a "hgnc_id" column, which contains a HGNC ID of a given gene. We then ask
# BioMart to return all mappings from that column to the `query_column`. For example, it can be the "ensembl_gene_id"
# column, which contains stable Ensembl gene ID.
BIOMART_REQUEST_TEMPLATE = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "{key_column}" value = "{identifier_list}"/>
        <Attribute name = "{key_column}" />
        <Attribute name = "{query_column}" />
    </Dataset>
</Query>""".replace('\n', '')

# This is to avoid getting a "URL too long" exception.
MAX_REQUEST_LENGTH = 5000


# Since we are using an external API call here, the @retry decorator will ensure that any sporadic network errors will
# be handled and the request will be retried.
@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def process_biomart_request(query):
    result = requests.get(query)
    # If there was an HTTP error, raise an exception. This will be caught by @retry.
    result.raise_for_status()
    return result.text


def split_into_chunks(lst, max_size, delimiter=','):
    """Splits an input list into chunks so that the total str representation of every chunk (including the delimiters)
    is no more than max_size."""
    lst = lst[:]
    chunks = []
    while lst:
        chunk_str_len = len(lst[0])
        chunk = [lst.pop(0)]
        while lst:
            size_incr = len(delimiter) + len(lst[0])
            if chunk_str_len + size_incr <= max_size:
                chunk_str_len += size_incr
                chunk.append(lst.pop(0))
            else:
                break
        chunks.append(chunk)
    return chunks


def query_biomart(key_column, query_column, identifier_list):
    """Query Ensembl BioMart with a list of identifiers (`identifier_list`) from one column (`key_column`) and return
    all mappings from those identifiers to another column (`query_column`) in form of a two-column Pandas dataframe.

    Args:
        key_column: A tuple of key column names in Ensembl and in the resulting dataframe, e.g. ('hgnc_id', 'HGNC_ID')
        query_column: A tuple of query column names, similar to `key_column`, e.g. ('ensembl_gene_id', 'EnsemblGeneID')
        identifier_list: List of identifiers to query, e.g. ['HGNC:10548', 'HGNC:10560']

    Returns:
        A Pandas dataframe with two columns. It will contain at most one row per input identifier. The query column will
        always contain a *list* to support the possibility of multiple mappings. In the example above, this will be:
               HGNC_ID      EnsemblGeneID
            0  HGNC:10548   [ENSG00000124788]
            1  HGNC:10560   [ENSG00000285258, ENSG00000163635]"""
    biomart_key_column, df_key_column = key_column
    biomart_query_column, df_query_column = query_column
    result = ''
    for identifier_chunk in split_into_chunks(identifier_list, MAX_REQUEST_LENGTH-len(BIOMART_REQUEST_TEMPLATE)):
        logger.info(f'Processing chunk of {len(identifier_chunk)} records')
        # Construct BioMart query from the template (see explanation above).
        biomart_query = BIOMART_REQUEST_TEMPLATE.format(
            key_column=biomart_key_column,
            query_column=biomart_query_column,
            identifier_list=','.join(identifier_chunk)
        )
        result += process_biomart_request(biomart_query)
    resulting_df = pd.read_table(StringIO(result), names=(df_key_column, df_query_column))
    # Group all potential mappings into lists.
    resulting_df = resulting_df.groupby(df_key_column)[df_query_column].apply(list).reset_index(name=df_query_column)
    return resulting_df
