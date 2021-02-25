#!/bin/bash
cd $1
input_dir_1=$2
input_dir_2=$3
result_dir=$4
true_gene_threshold=$5
Rscript Visualize_real_simulation_meta_analysis.R ${input_dir_1} ${input_dir_2} ${result_dir} ${true_gene_threshold}
rm -rf ~/.local/share/Trash/*

