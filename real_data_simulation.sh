#!/bin/bash
cd $1
input_dir=$2
result_dir=$3
true_gene_threshold=$4


for k in {1..15}
do
	echo $k
	Rscript real_data_simulation.R ${input_dir} ${result_dir} ${true_gene_threshold} ${k}
   	rm -rf ~/.local/share/Trash/*
done
