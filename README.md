# compareMetaDE
Integrated metaDE analysis with metapro and customized metaDE(v.1.0.5) package

## customized MetaDE package
Ordmeta and wFisher methods from metapro R package are integrated to MetaDE package to compare with.
You can download metapro package here. Download: <a href="https://github.com/unistbig/metapro/">https://github.com/unistbig/metapro/</a><br>
After metapro installation, you can install the the <i>metapro</i> package by typing 
open R and install the <i>customized_MetaDE</i> package by typing 

```
install.packages(file.path(path_to_file,"MetaDE_1.0.5_customized_bukyung.tar.gz"), repos = NULL, type="source")
```


## simulation data generation
run simulation_data_generation.sh on shell.

```
./simulation_data_generation.sh 'working directory where rnaseq_voom_norm_simulation.R exists' 'output_directory'
```
Count and voom folders will be created after the run, data in the voom folder will be used for meta-analysis.

