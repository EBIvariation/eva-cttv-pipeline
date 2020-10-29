#!/usr/bin/env python3

import gzip
from io import BytesIO
import os
import pandas
import pandas.testing
import pytest
import tempfile

from repeat_expansion_variants import pipeline


def get_test_resource(resource_name):
    """Gets full path to the test resource located in the same directory as the test module."""

    # E. g. repeat_expansion_variants/test/integration_test.py::test_pipeline (call)
    pytest_current_test = os.getenv('PYTEST_CURRENT_TEST')

    # E.g. repeat_expansion_variants/test/integration_test.py
    test_module_path = pytest_current_test.split(':')[0]

    # E. g. repeat_expansion_variants/test/
    test_module_dir = os.path.split(test_module_path)[0]

    print(os.path.join(test_module_dir, resource_name))

    # E. g.  repeat_expansion_variants/test/input_variant_summary.tsv
    return os.path.join(test_module_dir, resource_name)


def load_annotated_dataframe(dataframe_path, column_names=None):
    """
    Loads the dataframe from a TSV dump and prepares it for comparison:
        * Removes annotations (comments) and empty lines
        * Sorts rows by the values of all columns
        * Resets the index

    See example of an annotated dataframe in input_variant_summary.tsv.
    """

    # Load, remove annotations and blank lines
    df = pandas.read_table(dataframe_path, delimiter='\t', comment='#', skip_blank_lines=True, names=column_names)

    # Sort all rows using all columns
    df = df.sort_values(by=list(df.columns))

    # Reset index
    df = df.reset_index(drop=True)

    return df


@pytest.mark.skip(reason='Repeat expansion pipeline is not fully refactored yet')
def test_pipeline():
    """
    This is the integration test for the repeat expansion pipeline. It includes running the pipeline on a test dataset.
    The dataset (input_variant_summary.tsv) includes all ClinVar's "NT expansion" records as of 2020-04-08. The records
    of special interest are annotated and listed on top of that file. They include examples of all data peculiarities
    known to date:
        * All types of identifiers: HGVS-like coding/genomic; HGVS-like protein; Free text; non-standard; empty
        * Some variants having two records for coordinates in different genomic builds
        * Multiple associations with RCV; with GeneSymbol; and with both
        * Presence of different fields to determine repeat length: only sequence; only span; or both
        * Different repeat length: trinucleotide and non-trinucleotide
        * Repeat unit sequence including a N

    After running the pipeline, files produced by it are compared to the expected ones, which were manually examined
    for correctness.
    """

    # Load test dataset and remove annotations
    annotated_input_data = get_test_resource('input_variant_summary.tsv')
    input_data = ''.join([
         line
         for line in open(annotated_input_data)
         if line.strip() and not line.startswith('#')
    ])

    # Prepare input data, output paths, and run the pipeline
    input_file = BytesIO(gzip.compress(input_data.encode()))
    output_consequences, output_dataframe = [tempfile.NamedTemporaryFile(delete=False).name for _ in range(2)]
    pipeline.main(
        clinvar_xml=input_file,
        output_consequences=output_consequences,
        output_dataframe=output_dataframe
    )

    # Compare the dataframe outputs
    df_expected = load_annotated_dataframe(get_test_resource('output_dataframe.tsv'))
    df_actual = load_annotated_dataframe(output_dataframe)
    pandas.testing.assert_frame_equal(df_expected, df_actual)

    # Compare the consequence outputs. Since the actual output files do not contain headers, we specify the names
    # explicitly here
    consequences_column_names = ('RCVaccession', 'PlaceholderOnes', 'EnsemblGeneID',
                                 'EnsemblGeneName', 'RepeatType', 'PlaceholderZeroes')
    consequences_expected = load_annotated_dataframe(get_test_resource('output_consequences.tsv'),
                                                     consequences_column_names)
    consequences_actual = load_annotated_dataframe(output_consequences, consequences_column_names)
    pandas.testing.assert_frame_equal(consequences_expected, consequences_actual)
