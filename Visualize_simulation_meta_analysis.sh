#!/bin/bash
total_study=20
i=(1)
cd $1
input_dir=$2
result_dir=$3

Rscript Visualize_simulation_meta_analysis.R ${i} ${total_study} ${input_dir} ${result_dir}
rm -rf ~/.local/share/Trash/*

