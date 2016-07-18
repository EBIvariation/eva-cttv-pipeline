import json
import gzip
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
import http.client
from collections import UserDict

from eva_cttv_pipeline import utilities


__author__ = 'Javier Lopez: javild@gmail.com'


def get_rcv_to_rsnsv_mapping(variant_summary_file):

    rcv_to_rs = {}
    rcv_to_nsv = {}

    print('Loading mapping RCV->rs/nsv')
    fdr = utilities.open_file(variant_summary_file, "rt")
    fdr.readline()
    for line in fdr:
        parts = line.split('\t')
        rcv_list = parts[8].split(';')
        if parts[6] != '-' and parts[6] != '-1' and parts[6] != '':
            for rcv_id in rcv_list:
                rcv_to_rs[rcv_id] = 'rs' + parts[6]
        if parts[7] != '-' and parts[7] != '-1' and parts[7] != '':
            for rcv_id in rcv_list:
                rcv_to_nsv[rcv_id] = parts[7]
    fdr.close()
    print(' Done.')

    return rcv_to_rs, rcv_to_nsv


class ClinvarRecord(UserDict):
    """
    Class of which instances hold data on individual clinvar records. Subclass of UserDict rather
    than dict in order to use attributes
    """

    cached_symbol_2_ensembl = {}

    score_map = {
        "CLASSIFIED_BY_SINGLE_SUBMITTER": 1,
        "NOT_CLASSIFIED_BY_SUBMITTER": None,
        "CLASSIFIED_BY_MULTIPLE_SUBMITTERS": 2,
        "REVIEWED_BY_EXPERT_PANEL": 3,
        "REVIEWED_BY_PROFESSIONAL_SOCIETY": 4
    }

    def __init__(self, mappings=None, a_dictionary=None):
        UserDict.__init__(self, dict=a_dictionary)
        self.rs = self.__get_rs(mappings.rcv_to_rs)
        self.nsv = self.__get_nsv(mappings.rcv_to_nsv)
        self.consequence_type = self.__get_main_consequence_types(mappings.consequence_type_dict,
                                                                  mappings.rcv_to_rs)

    @property
    def gene_id(self):
        j = 0
        measure = self.data['referenceClinVarAssertion']['measureSet']['measure']
        found = False
        while j < len(measure) and not found:
            attribute_set = measure[j]['attributeSet']
            i = 0
            while i < len(attribute_set) and \
                    not attribute_set[i]['attribute']['type'].startswith('HGVS'):
                i += 1
            found = (i < len(attribute_set))
            j += 1

        if found:
            return attribute_set[i]['attribute']['value'][:9]
        else:
            return self.data['referenceClinVarAssertion']['measureSet']['measure'][0]['name'][0]['elementValue']['value']

    @property
    def ensembl_id(self):
        measure = self.data['referenceClinVarAssertion']['measureSet']['measure']
        for j in range(len(measure)):
            if 'measureRelationship' not in measure[j]:
                continue
            measure_relationship = measure[j]['measureRelationship']
            for l in range(len(measure_relationship)):
                if 'symbol' not in measure_relationship[l]:
                    continue
                symbol = measure_relationship[l]['symbol']
                for i in range(len(symbol)):
                    try:
                        return ClinvarRecord.cached_symbol_2_ensembl[
                            symbol[i]['elementValue']['value']
                        ]
                    except KeyError:
                        not_solved = True
                        while not_solved:
                            try:
                                url = 'http://rest.ensembl.org/' + \
                                      'lookup/symbol/homo_sapiens/' + \
                                      symbol[i]['elementValue']['value'] + \
                                      '?content-type=application/json'
                                raw_reply = urllib.request.urlopen(url).read()
                                ensembl_json = json.loads(raw_reply.decode())
                                not_solved = False
                            except urllib.error.HTTPError as e:
                                if e.code == 400:
                                    print('WARNING: Bad request code returned +' +
                                          'from ENSEMBL rest.')
                                    print(' ClinVar accession: ' + self.accession)
                                    print(' Gene symbol: ' +
                                          symbol[i]['elementValue']['value'])
                                    print(' Error: ')
                                    print(e)
                                    print(' Returning None.')
                                    ClinvarRecord.cached_symbol_2_ensembl[
                                        symbol[i]['elementValue']['value']] = None
                                    return None
                                else:
                                    time.sleep(0.05)
                            except http.client.BadStatusLine:
                                time.sleep(3)

                        if len(ensembl_json) > 0:
                            if ensembl_json['id'] and \
                                    (ensembl_json['object_type'] == 'Gene'):
                                ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = \
                                    ensembl_json['id']
                                return ensembl_json['id']
                            else:
                                print("WARNING at clinvar_record.py: ENSEMBL's REST API " +
                                      "returned an unexpected json object")
                                print(" Queried gene symbol: " +
                                      symbol[i]['elementValue']['value'])
                                print(" ENSEMBL's API response:")
                                print(ensembl_json)
                                print("Exiting.")
                                sys.exit(1)
                        else:
                            ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = \
                                None
        return None

    @property
    def date(self):
        return datetime.fromtimestamp(
            self.data['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    @property
    def score(self):
        return self.score_map[
            self.data['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']
        ]

    @property
    def accession(self):
        return self.data['referenceClinVarAssertion']['clinVarAccession']['acc']

    @property
    def traits(self):
        trait_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            trait_list.append([])
            for name in trait['name']:
                # First trait name in the list will always be the "Preferred" one
                if name['elementValue']['type'] == 'Preferred':
                    trait_list[-1] = [name['elementValue']['value']] + trait_list[-1]
                elif name['elementValue']['type'] in ["EFO URL", "EFO id", "EFO name"]:
                    continue  # if the trait name not originally from clinvar
                else:
                    trait_list[-1].append(name['elementValue']['value'])

        return trait_list

    @property
    def trait_pubmed_refs(self):
        pubmed_refs_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            pubmed_refs_list.append([])
            if 'citation' in trait:
                for citation in trait['citation']:
                    if ('id' in citation) and citation['id'] is not None:
                        if citation['id']['source'] == 'PubMed':
                            pubmed_refs_list[-1].append(int(citation['id']['value']))

        return pubmed_refs_list

    @property
    def observed_pubmed_refs(self):
        pubmed_refs_list = []
        if 'observedIn' in self.data['referenceClinVarAssertion']:
            for observed_in in self.data['referenceClinVarAssertion']['observedIn']:
                for observed_data in observed_in['observedData']:
                    if 'citation' in observed_data:
                        for citation in observed_data['citation']:
                            if ('id' in citation) and citation['id'] is not None:
                                if citation['id']['source'] == 'PubMed':
                                    pubmed_refs_list.append(int(citation['id']['value']))
        return pubmed_refs_list

    @property
    def measure_set_pubmed_refs(self):
        pubmed_refs_list = []
        for measure in self.data['referenceClinVarAssertion']['measureSet']['measure']:
            if 'citation' in measure:
                for citation in measure['citation']:
                    if 'id' in citation and citation['id'] is not None:
                        if citation['id']['source'] == 'PubMed':
                            pubmed_refs_list.append(int(citation['id']['value']))
        return pubmed_refs_list

    @property
    def trait_refs_list(self):
        return [['http://europepmc.org/abstract/MED/' + str(ref) for ref in ref_list]
                for ref_list in self.trait_pubmed_refs]

    @property
    def observed_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.observed_pubmed_refs]

    @property
    def measure_set_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.measure_set_pubmed_refs]

    @property
    def hgvs(self):
        hgvs_list = []
        for measure in self.data['referenceClinVarAssertion']['measureSet']['measure']:
            for attribute_set in measure['attributeSet']:
                if attribute_set['attribute']['type'].startswith('HGVS'):
                    hgvs_list.append(attribute_set['attribute']['value'])

        return hgvs_list

    @property
    def clinical_significance(self):
        return \
            self.data['referenceClinVarAssertion']['clinicalSignificance']['description'].lower()

    def __get_rs(self, rcv_to_rs):
        try:
            return rcv_to_rs[self.accession]
        except KeyError:
            return None

    def __get_nsv(self, rcv_to_nsv):
        try:
            return rcv_to_nsv[self.accession]
        except KeyError:
            return None

    def __get_main_consequence_types(self, consequence_type_dict, rcv_to_rs):

        new_rs_id = self.__get_rs(rcv_to_rs)

        if new_rs_id is not None and new_rs_id in consequence_type_dict:
            return consequence_type_dict[new_rs_id]
        elif self.accession in consequence_type_dict:
            return consequence_type_dict[self.accession]
        else:
            return None

    @property
    def variant_type(self):
        return self.data['referenceClinVarAssertion']['measureSet']['measure'][0]['type']

    @property
    def allele_origins(self):
        allele_origins = set()
        for clinvar_assetion_document in self.data['clinVarAssertion']:
            for observed_in_document in clinvar_assetion_document['observedIn']:
                allele_origins.add(observed_in_document['sample']['origin'])

        return list(allele_origins)
