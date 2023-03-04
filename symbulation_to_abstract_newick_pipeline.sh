#! /usr/bin/bash

python abstract_phylogenies.py $1 $2 $3
python make_links.py abstract_interactions.csv

symfile=`basename $1`
hostfile=`basename $2`

alifedata-phyloinformatics-convert fromalifedata --input-file compressed_sym.csv --output-file symtree.nwk --output-schema newick --suppress-unifurcations
alifedata-phyloinformatics-convert fromalifedata --input-file compressed_host.csv --output-file hosttree.nwk --output-schema newick --suppress-unifurcations