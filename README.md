# compareMetaDE
Integrative pipeline with metapro, customized metaDE(v.1.0.5) package and data generation steps for DE gene meta-analysis.

## customized MetaDE package
Ordmeta and wFisher methods from metapro R package are integrated to MetaDE package to compare with.
You can download metapro package here. Download: <a href="https://github.com/unistbig/metapro/">https://github.com/unistbig/metapro/</a><br>
After metapro installation, you can install the the <i>metapro</i> package by typing 
open R and install the <i>customized_MetaDE</i> package by typing 

```
install.packages(file.path(path_to_file,"MetaDE_1.0.5_customized_bukyung.tar.gz"), repos = NULL, type="source")
```
The original MetaDE package(v.1.0.5) is avaialable at CRAN. 
If you want to check, Download: <a href=" https://cran.r-project.org/src/contrib/Archive/MetaDE/MetaDE_1.0.5.tar.gz/"> https://cran.r-project.org/src/contrib/Archive/MetaDE//</a><br>

## :paperclip: Dependency
check Dependencies before running main analysis.

* [tidyverse](https://github.com/tidyverse/tidyverse)
* [magrittr](https://github.com/tidyverse/magrittr)

### RNA-seq simulation
* [Biobase](https://bioconductor.org/packages/Biobase)
* [edgeR](https://bioconductor.org/packages/edgeR)
* [limma](https://bioconductor.org/packages/limma)
* [SimSeq](https://CRAN.R-project.org/package=SimSeq)

### meta analysis
* [metapro](https://github.com/unistbig/metapro)
* [MetaQC](https://CRAN.R-project.org/package=MetaQC) -install archived version, 0.1.13
* [MetaDE](https://github.com/unistbig/compareMetaDE/blob/main/MetaDE_1.0.5_customized_bukyung.tar.gz) -install customized package I uploaded
 
 ### Visualization
* [ggplot2](https://CRAN.R-project.org/package=ggplot2)
* [reshape2](https://cran.r-project.org/package=reshape2)
* [ROCR](https://cran.r-project.org/package=ROCR)
* [colortools](https://cran.r-project.org/package=colortools)
* [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)
* [VennDiagram](https://CRAN.R-project.org/package=VennDiagram)

## Simualation data analysis
### simulation data generation
run simulation_data_generation.sh on shell.

```
./simulation_data_generation.sh 'working directory where rnaseq_voom_norm_simulation.R exists' 'output_directory'
```
Count and voom folders will be created after the run, data in the voom folder will be used for meta-analysis.

### simulation data meta-analysis
run run_meta_simul.sh on shell.
```
./simulation_data_generation.sh 'working directory where Simulation_meta_analysis.R exists' 'input_directory where voom folder locate' 'output_directory'
```
 Results will be saved in "output_directory/Meta_Res" folder.

### Visualize simulation meta-analysis
run Visualize_simulation_meta_analysis.sh on shell.
```
./Visualize_simulation_meta_analysis.sh 'working directory where Visualize_simulation_meta_analysis.R exists' 'input_directory where voom folder locate' 'output_directory'
```
 Plots will be saved in "output_directory/Simulation_analysis_plot" folder.
 
 ## Microarray data analysis
 ### True gene selection
 run true_gene_selection.sh on shell.
```
./true_gene_selection.sh 'working directory where true_gene_selection.R exists' 'input_directory where prostate data locate' 'output_directory' 'true_de_gene_threshold'
```
 Truegenes will be saved in output_directory.
 Other 6 study and 9 study analysis results will be saved in "output_directory/Meta_Res_real" folder.
 We used 0.01 of FDR threshold for true de gene selection.
 
 ### real data simulation
 run real_data_simulation.sh on shell
 ```
./real_data_simulation.sh 'working directory where real_data_simulation.R exists' 'input_directory where prostate data locate' 'output_directory' 'true_de_gene_threshold'
```
 Results will be saved in "output_directory/Meta_Res_real" folder.
 
 ### Visualize real data analysis
 run Visualize_real_simulation_meta_analysis.sh on shell
 ```
./Visualize_real_simulation_meta_analysis.sh 'working directory where Visualize_real_simulation_meta_analysis.R exists' 'input_directory1 where real data simulation locate(Meta_Res_real folder we made)' 'input_directory2 where selected true genes locate' 'output_directory' 'true_de_gene_threshold'
```
 Results will be saved in "output_directory/Real_study_analysis_plot" folder.
 
 
