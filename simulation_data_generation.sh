#!/bin/bash
n_diffexps=(100 300 600)
total_study=20
associated_study_nums=(2 5 10)
is=(1 2 3 4 5 6 7 8 9 10)
cd $1
result_dir=$2
for i in ${is[@]}
do
   for n_diffexp in ${n_diffexps[@]}
   do
   	for associated_study_num in ${associated_study_nums[@]}
	do
		echo $i
		echo $n_diffexp
		echo $associated_study_num

		Rscript rnaseq_voom_norm_simulation.R ${i} ${n_diffexp} ${total_study} ${associated_study_num} ${result_dir}
   		rm -rf ~/.local/share/Trash/*
   	done
   done
done

