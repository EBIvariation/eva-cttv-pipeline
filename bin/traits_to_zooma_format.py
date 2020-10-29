#!/usr/bin/env python3

import argparse
import itertools
from time import gmtime, strftime

from eva_cttv_pipeline import clinvar_xml_utils


class OntologyUri:
    db_to_uri_dict = {
        'orphanet': 'http://www.orpha.net/ORDO/Orphanet_{}',
        'omim':     'http://identifiers.org/omim/{}',
        'efo':      'http://www.ebi.ac.uk/efo/{}',
        'mesh':     'http://identifiers.org/mesh/{}',
        'medgen':   'http://identifiers.org/medgen/{}',
        'mondo':    'http://purl.obolibrary.org/obo/MONDO_{}',
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri


def write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, date, outfile):
    zooma_output_list = [clinvar_acc,
                         variant_id,
                         'disease',
                         trait_name,
                         str(ontology_uri),
                         'clinvar-xrefs',
                         date]
    outfile.write('\t'.join(zooma_output_list) + '\n')


def process_clinvar_record(clinvar_record, outfile):
    if clinvar_record.measure is None:
        return
    variant_ids = [variant_id
                   for variant_id in (clinvar_record.measure.rs_id, clinvar_record.measure.nsv_id)
                   if variant_id is not None]
    traits = clinvar_record.traits
    for variant_id, trait in itertools.product(variant_ids, traits):
        if trait.name.lower() == 'not provided':
            continue
        for db, identifier, status in trait.xrefs:
            if status != 'current' or db.lower() not in OntologyUri.db_to_uri_dict:
                continue
            ontology_uri = OntologyUri(identifier, db)
            write_zooma_record(clinvar_record.accession, variant_id, trait.name, ontology_uri,
                               strftime('%d/%m/%y %H:%M', gmtime()), outfile)


def main(clinvar_xml, zooma_feedback):
    with open(zooma_feedback, 'wt') as outfile:
        outfile.write('STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE\n')
        for clinvar_record in clinvar_xml_utils.ClinVarDataset(clinvar_xml):
            process_clinvar_record(clinvar_record, outfile)


parser = argparse.ArgumentParser('Extracts OMIM and MedGen cross-references from ClinVar for submission to ZOOMA')
parser.add_argument('--clinvar-xml',    help='ClinVar XML release',                 required=True)
parser.add_argument('--zooma-feedback', help='Disease string to ontology mappings', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.clinvar_xml, args.zooma_feedback)
