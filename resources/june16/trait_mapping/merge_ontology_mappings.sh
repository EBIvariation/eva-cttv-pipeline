#!/usr/bin/env bash

# trait_to_urls file contains only traits that have mappings
# feb_2016_efo_maps.csv is file containing mappings from feb16 batch

trait_to_urls_file=$1
old_mappings_file=$2
output_file=$3

grep -v -f <(cut -f1 $trait_to_urls_file | sed -e 's/$/\t/') $old_mappings_file | grep http > $output_file
cat $trait_to_urls_file >> $output_file
