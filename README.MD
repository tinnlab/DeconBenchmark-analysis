This repository contains code to reproduce the benchmarking results in the article "Twelve Years of Cellular Deconvolution: Applications, Benchmark, Methodology, and Challenges".

For standalone R package to run deconvolution methods, visit http://deconvolution.us.to/gitlab/deconvolution/benchmark.

To get started, follow the steps below:

1. Clone this repository
2. Go to https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5, 
download the `rds` file of `Tabula Sapiens - All Cells`, rename it to `local.rds` and place it in the `data` directory.
3. Install docker:
   Follow instructions at https://docs.docker.com/engine/install/ to install docker.
   All prebuilt docker images can be found at https://hub.docker.com/u/deconvolution.
   Docker will automatically pull the images when running the scripts below.
   However, CIBERSORT has a specific license requirement and is not available on docker hub.
   Please visit http://deconvolution.us.to/gitlab/deconvolution/docker to read more about how to build the CIBERSORT docker image
   and any other method.
4. Using `R` version >= 4.0.0, run the following commands to install necessary packages used in the analysis:
```r
if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install(c("tidyverse", "babelwhale", "RcppHungarian", "stringr", "processx", "rhdf5", "Seurat", "edgeR"))
```
5. Visit the Matlab license center at http://www.mathworks.com/licensecenter/licenses to obtain a license file.
The hostid of this license is the same as the hostid of the host machine. The user of this license must be root.
This license is required only when running BayesCCE, Deblender, and DecOT.
Save the license file and replace `PATH_TO_MATLAB_LICENSE` in R scripts below with the ABSOLUTE path to the license file.
6. Run `generate-simulation-data.R` to generate simulated data.
All simulated data will be saved in:
   - `data/sim` for accuracy and consistency benchmarking
   - `data/scalability` for scalability benchmarking
7. Run `accuracy.R` to generate the accuracy results.
All results will be saved in the `results/accuracy` directory, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>.rds`.
   Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`. 
   This `rds` file contains a dataframe with method name, dataset name, and 4 accuracy metrics.
8. Run `consistency.R` to generate the consistency results.
   All results will be saved in the `results/consistency` directory, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>_<seed>.rds`.
     Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`.
     This `rds` file contains a dataframe with method name, dataset name, sample, seed, MAE, and correlation.
9. Run `scalability.R` to generate the scalability results.
At the same time, run `memory-monitor.R` to monitor the memory usage of the running docker containers.
All scalability results will be saved in the `results/scalability` directory, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>_<numberOfBulkSample>.rds`.
     Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`.
     This `rds` file contains a dataframe with method name, dataset name, number of bulk samples, running time, and memory usage.


#### Session Info
```
R version 4.1.3 (2022-03-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS/LAPACK: /home/..../lib/libopenblasp-r0.3.21.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RcppHungarian_0.2  rhdf5_2.38.1       babelwhale_1.1.0   SeuratObject_4.1.3 Seurat_4.3.0       lubridate_1.9.2   
 [7] forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0        purrr_1.0.1        readr_2.1.4        tidyr_1.3.0       
[13] tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0   

loaded via a namespace (and not attached):
  [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
  [6] rprojroot_2.0.3        spatstat.data_3.0-0    leiden_0.4.3           listenv_0.9.0          remotes_2.4.2         
 [11] dynutils_1.0.11        ggrepel_0.9.3          fansi_1.0.4            codetools_0.2-19       splines_4.1.3         
 [16] polyclip_1.10-4        jsonlite_1.8.4         ica_1.0-3              cluster_2.1.4          png_0.1-8             
 [21] uwot_0.1.14            shiny_1.7.4            sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.1.3        
 [26] httr_1.4.5             assertthat_0.2.1       Matrix_1.5-3           fastmap_1.1.1          lazyeval_0.2.2        
 [31] limma_3.50.3           cli_3.6.0              later_1.3.0            htmltools_0.5.4        tools_4.1.3           
 [36] igraph_1.4.1           gtable_0.3.1           glue_1.6.2             RANN_2.6.1             reshape2_1.4.4        
 [41] Rcpp_1.0.10            scattermore_0.8        rhdf5filters_1.6.0     vctrs_0.5.2            spatstat.explore_3.0-6
 [46] nlme_3.1-162           progressr_0.13.0       lmtest_0.9-40          spatstat.random_3.1-3  ps_1.7.2              
 [51] globals_0.16.2         timechange_0.2.0       mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
 [56] irlba_2.3.5.1          goftest_1.2-3          future_1.31.0          edgeR_3.36.0           MASS_7.3-58.2         
 [61] zoo_1.8-11             scales_1.2.1           hms_1.1.2              promises_1.2.0.1       spatstat.utils_3.0-1  
 [66] RColorBrewer_1.1-3     reticulate_1.28        pbapply_1.7-0          gridExtra_2.3          stringi_1.7.12        
 [71] desc_1.4.2             rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0     lattice_0.20-45       
 [76] Rhdf5lib_1.16.0        ROCR_1.0-11            tensor_1.5             patchwork_1.1.2        htmlwidgets_1.6.1     
 [81] processx_3.8.0         cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.34.0      RcppAnnoy_0.0.20      
 [86] plyr_1.8.8             magrittr_2.0.3         R6_2.5.1               generics_0.1.3         pillar_1.8.1          
 [91] withr_2.5.0            proxyC_0.3.3           fitdistrplus_1.1-8     survival_3.5-3         abind_1.4-5           
 [96] sp_1.6-0               future.apply_1.10.0    crayon_1.5.2           KernSmooth_2.23-20     utf8_1.2.3            
[101] spatstat.geom_3.0-6    plotly_4.10.1          tzdb_0.3.0             locfit_1.5-9.7         grid_4.1.3            
[106] data.table_1.14.8      digest_0.6.31          xtable_1.8-4           httpuv_1.6.9           RcppParallel_5.1.6    
[111] munsell_0.5.0          viridisLite_0.4.1     
```
