from collections import defaultdict
import copy


def main():
    rcv_to_urls = get_rcv_to_urls("/home/tom/Job_Working_Directory/Python/eva_cttv_pipeline/resources/june16/trait_mapping/eva_clinvar_mappings_collated_july_2016_sheet2_col2-3.csv")

    rcv_to_trait = get_rcv_to_trait("/home/tom/Job_Working_Directory/Python/eva_cttv_pipeline/resources/june16/trait_mapping/eva_clinvar_mappings_collated_july_2016_sheet1_col1-2.csv")

    trait_to_urls, rcvs_no_url, rcvs_no_trait, trait_no_url = get_trait_to_url(rcv_to_trait, rcv_to_urls)

    with open("/home/tom/Job_Working_Directory/Python/eva_cttv_pipeline/resources/june16/trait_mapping/trait_to_urls.tsv", "wt") as f:
        for trait, urls in trait_to_urls.items():
            url_string = ','.join(urls)
            out_string = '\t'.join([trait, url_string])
            f.write(out_string + "\n")

    print(len(trait_no_url))
    print(len(trait_to_urls))


def get_rcv_to_urls(file_path):
    rcv_to_urls = {}
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            rcv = line_list[0]
            urls_str = line_list[1]
            urls_set = set(urls_str.split(","))
            rcv_to_urls[rcv] = urls_set
    return rcv_to_urls


def get_rcv_to_trait(file_path):
    rcv_to_trait = {}
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            rcv = line_list[0]
            trait = line_list[1]
            rcv_to_trait[rcv] = trait
    return rcv_to_trait


def get_trait_to_url(rcv_to_trait, rcv_to_urls):
    rcv_to_trait = copy.deepcopy(rcv_to_trait)
    rcv_to_urls = copy.deepcopy(rcv_to_urls)

    rcvs_no_url = set()
    rcvs_no_trait = set()
    trait_no_url = set()
    trait_to_urls = defaultdict(set)
    for rcv, trait in copy.deepcopy(rcv_to_trait).items():
        if rcv not in rcv_to_urls:
            rcvs_no_url.add(rcv)
            trait_no_url.add(trait)
            del rcv_to_trait[rcv]
            continue

        new_urls = rcv_to_urls[rcv]
        trait_to_urls[trait].update(new_urls)

        del rcv_to_urls[rcv]

        if trait in trait_no_url:
            trait_no_url.remove(trait)

    for rcv, urls in rcv_to_urls.items():
        rcvs_no_trait.add(rcv)

    return trait_to_urls, rcvs_no_url, rcvs_no_trait, trait_no_url


if __name__ == "__main__":
    main()
