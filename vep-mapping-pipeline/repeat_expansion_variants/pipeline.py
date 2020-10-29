#!/usr/bin/env python3
"""A pipeline to extract repeat expansion variants from ClinVar TSV dump. For documentation refer to README.md"""

import logging

import numpy as np
import pandas as pd

from eva_cttv_pipeline import clinvar_xml_utils
from . import biomart, clinvar_identifier_parsing

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


STANDARD_CHROMOSOME_NAMES = {str(c) for c in range(1, 23)} | {'X', 'Y', 'M', 'MT'}


def none_to_nan(*args):
    """Converts all arguments which are None to np.nan, for consistency inside a Pandas dataframe."""
    return [np.nan if a is None else a for a in args]


def load_clinvar_data(clinvar_xml):
    """Load ClinVar data, preprocess, and return it as a Pandas dataframe."""
    # Iterate through ClinVar XML records. Load all Microsatellite variants
    variant_data = []  # Name, RCVaccession, GeneSymbol, HGNC_ID
    for i, clinvar_record in enumerate(clinvar_xml_utils.ClinVarDataset(clinvar_xml)):
        if i % 10000 == 0:
            logger.info(f'Processed {i} records, have {len(variant_data)} variants')
        if clinvar_record.measure and clinvar_record.measure.variant_type == 'Microsatellite':
            # Extract gene symbol(s). Here and below, dashes are sometimes assigned to be compatible with the variant
            # summary format which was used previously.
            gene_symbols = clinvar_record.measure.preferred_gene_symbols
            if not gene_symbols:
                gene_symbols = ['-']
            # Extract HGNC ID
            hgnc_ids = clinvar_record.measure.hgnc_ids
            hgnc_id = hgnc_ids[0] if len(hgnc_ids) == 1 and len(gene_symbols) == 1 else '-'
            # Append data strings
            for gene_symbol in gene_symbols:
                variant_data.append([
                    clinvar_record.measure.name,
                    clinvar_record.accession,
                    gene_symbol,
                    hgnc_id
                ])
    variants = pd.DataFrame(variant_data, columns=('Name', 'RCVaccession', 'GeneSymbol', 'HGNC_ID'))
    # Since the same record can have coordinates in multiple builds, it can be repeated. Remove duplicates
    variants = variants.drop_duplicates()
    # Sort values by variant name
    return variants.sort_values(by=['Name'])


def parse_variant_identifier(row):
    """Parse variant identifier and extract certain characteristics into separate columns."""
    variant_name = str(row.Name)
    row['TranscriptID'], row['CoordinateSpan'], row['RepeatUnitLength'], row['IsProteinHGVS'] = \
        none_to_nan(*clinvar_identifier_parsing.parse_variant_identifier(variant_name))
    return row


def annotate_ensembl_gene_info(variants):
    """Annotate the `variants` dataframe with information about Ensembl gene ID and name"""

    # Ensembl gene ID can be determined using three ways, listed in the order of decreasing priority. Having multiple
    # ways is necessary because no single method works on all ClinVar variants.
    gene_annotation_sources = (
        # Dataframe column    Biomart column   Filtering function
        ('HGNC_ID',                 'hgnc_id', lambda i: i.startswith('HGNC:')),
        ('GeneSymbol',   'external_gene_name', lambda i: i != '-'),
        ('TranscriptID',        'refseq_mrna', lambda i: pd.notnull(i)),
    )
    # This copy of the dataframe is required to facilitate filling in data using the `combine_first()` method. This
    # allows us to apply priorities: e.g., if a gene ID was already populated using HGNC_ID, it will not be overwritten
    # by a gene ID determined using GeneSymbol.
    variants_original = variants.copy(deep=True)

    for column_name_in_dataframe, column_name_in_biomart, filtering_function in gene_annotation_sources:
        # Get all identifiers we want to query BioMart with
        identifiers_to_query = sorted({
            i for i in variants[column_name_in_dataframe]
            if filtering_function(i)
        })
        # Query BioMart for Ensembl Gene IDs
        annotation_info = biomart.query_biomart(
            key_column=(column_name_in_biomart, column_name_in_dataframe),
            query_column=('ensembl_gene_id', 'EnsemblGeneID'),
            identifier_list=identifiers_to_query,
        )
        # Make note where the annotations came from
        annotation_info['GeneAnnotationSource'] = column_name_in_dataframe
        # Combine the information we received with the *original* dataframe (a copy made before any iterations of this
        # cycle were allowed to run). This is similar to SQL merge.
        annotation_df = pd.merge(variants_original, annotation_info, on=column_name_in_dataframe, how='left')
        # Update main dataframe with the new values. This replaces the NaN values in the dataframe with the ones
        # available in another dataframe we just created, `annotation_df`.
        variants = variants \
            .set_index([column_name_in_dataframe]) \
            .combine_first(annotation_df.set_index([column_name_in_dataframe]))

    # Reset index to default
    variants.reset_index(inplace=True)
    # Some records are being annotated to multiple Ensembl genes. For example, HGNC:10560 is being resolved to
    # ENSG00000285258 and ENSG00000163635. We need to explode dataframe by that column.
    variants = variants.explode('EnsemblGeneID')

    # Based on the Ensembl gene ID, annotate (1) gene name and (2) which chromosome it is on
    gene_query_columns = (
        ('external_gene_name', 'EnsemblGeneName'),
        ('chromosome_name', 'EnsemblChromosomeName'),
    )
    for column_name_in_biomart, column_name_in_dataframe in gene_query_columns:
        annotation_info = biomart.query_biomart(
            key_column=('ensembl_gene_id', 'EnsemblGeneID'),
            query_column=(column_name_in_biomart, column_name_in_dataframe),
            identifier_list=sorted({str(i) for i in variants['EnsemblGeneID'] if str(i).startswith('ENSG')}),
        )
        variants = pd.merge(variants, annotation_info, on='EnsemblGeneID', how='left')
        # Check that there are no multiple mappings for any given ID
        assert variants[column_name_in_dataframe].str.len().dropna().max() == 1, \
            'Found multiple gene ID → gene attribute mappings!'
        # Convert the one-item list into a plain column
        variants = variants.explode(column_name_in_dataframe)

    return variants


def determine_repeat_type(row):
    """Based on all available information about a variant, determine its type. The resulting type can be:
        * trinucleotide_repeat_expansion, corresponding to SO:0002165
        * short_tandem_repeat_expansion, corresponding to SO:0002162
        * NaN (not able to determine)
    Also, depending on the information, determine whether the record is complete, i.e., whether it has all necessary
    fields to be output for the final "consequences" table."""
    repeat_type = np.nan
    if row['IsProteinHGVS']:
        # For protein HGVS notation, assume that repeat is a trinucleotide one, since it affects entire amino acids
        repeat_type = 'trinucleotide_repeat_expansion'
    else:
        # As a priority, use the repeat unit length determined directly from base sequence
        repeat_unit_length = row['RepeatUnitLength']
        # If not available, fall back to using and end coordinate difference
        if pd.isnull(repeat_unit_length):
            repeat_unit_length = row['CoordinateSpan']
        # Determine repeat type based on repeat unit length
        if pd.notnull(repeat_unit_length):
            if repeat_unit_length % 3 == 0:
                repeat_type = 'trinucleotide_repeat_expansion'
            else:
                repeat_type = 'short_tandem_repeat_expansion'
    row['RepeatType'] = repeat_type
    # Based on the information which we have, determine whether the record is complete
    row['RecordIsComplete'] = (
        pd.notnull(row['EnsemblGeneID']) and
        pd.notnull(row['EnsemblGeneName']) and
        pd.notnull(row['RepeatType']) and
        row['EnsemblChromosomeName'] in STANDARD_CHROMOSOME_NAMES
    )
    return row


def generate_output_files(variants, output_consequences, output_dataframe):
    """Postprocess and output final tables."""

    # Rearrange order of dataframe columns
    variants = variants[
        ['Name', 'RCVaccession', 'GeneSymbol', 'HGNC_ID',
         'RepeatUnitLength', 'CoordinateSpan', 'IsProteinHGVS', 'TranscriptID',
         'EnsemblGeneID', 'EnsemblGeneName', 'EnsemblChromosomeName', 'GeneAnnotationSource',
         'RepeatType', 'RecordIsComplete']
    ]
    # Write the full dataframe. This is used for debugging and investigation purposes.
    variants.sort_values(by=['Name', 'RCVaccession', 'GeneSymbol'])
    variants.to_csv(output_dataframe, sep='\t', index=False)

    # Generate consequences table
    consequences = variants[variants['RecordIsComplete']] \
        .groupby(['RCVaccession', 'EnsemblGeneID', 'EnsemblGeneName'])['RepeatType'] \
        .apply(set).reset_index(name='RepeatType')
    # Check that for every (RCV, gene) pair there is only one consequence type
    assert consequences['RepeatType'].str.len().dropna().max() == 1, 'Multiple (RCV, gene) → variant type mappings!'
    # Get rid of sets
    consequences['RepeatType'] = consequences['RepeatType'].apply(list)
    consequences = consequences.explode('RepeatType')
    # Form a six-column file compatible with the consequence mapping pipeline, for example:
    # RCV000005966    1    ENSG00000156475    PPP2R2B    trinucleotide_repeat_expansion    0
    consequences['PlaceholderOnes'] = 1
    consequences['PlaceholderZeroes'] = 0
    consequences = consequences[['RCVaccession', 'PlaceholderOnes', 'EnsemblGeneID', 'EnsemblGeneName', 'RepeatType',
                                 'PlaceholderZeroes']]
    consequences.sort_values(by=['RepeatType', 'RCVaccession', 'EnsemblGeneID'], inplace=True)
    # Check that there are no empty cells in the final consequences table
    assert consequences.isnull().to_numpy().sum() == 0
    # Write the consequences table. This is used by the main evidence string generation pipeline.
    consequences.to_csv(output_consequences, sep='\t', index=False, header=False)


def main(clinvar_xml, output_consequences, output_dataframe):
    """Process data and generate output files.

    Args:
        clinvar_xml: filepath to the ClinVar XML file.
        output_consequences: filepath to the output file with variant consequences. The file uses a 6-column format
            compatible with the VEP mapping pipeline (see /vep-mapping-pipeline/README.md).
        output_dataframe: filepath to the output file with the full dataframe used in the analysis. This will contain
            all relevant columns and can be used for review or debugging purposes."""

    logger.info('Load and preprocess variant data')
    variants = load_clinvar_data(clinvar_xml)

    logger.info('Parse variant names and extract information about transcript ID and repeat length')
    variants = variants.apply(lambda row: parse_variant_identifier(row), axis=1)

    logger.info('Match each record to Ensembl gene ID and name')
    variants = annotate_ensembl_gene_info(variants)

    logger.info('Determine variant type and whether the record is complete')
    variants = variants.apply(lambda row: determine_repeat_type(row), axis=1)

    logger.info('Postprocess data and output the two final tables')
    generate_output_files(variants, output_consequences, output_dataframe)

    logger.info('Completed successfully')
