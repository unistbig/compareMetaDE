#!/bin/bash
cd $1
input_dir=$2
result_dir=$3
trun_gene_threshold=$4

Rscript true_gene_selection.R ${input_dir} ${result_dir} ${trun_gene_threshold}
rm -rf ~/.local/share/Trash/*
