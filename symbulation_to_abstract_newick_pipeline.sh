#! /usr/bin/bash

python abstract_phylogenies.py $1 $2 $3
python make_links.py abstract_interactions

symfile=`basename $1`
hostfile=`basename $2`

alifedata-phyloinformatics-convert fromalifedata --input-file compressed_${symfile} --output-file symtree.nwk --output-schema newick --suppress-unifurcations
alifedata-phyloinformatics-convert fromalifedata --input-file compressed_${hostfile} --output-file hosttree.nwk --output-schema newick --suppress-unifurcations