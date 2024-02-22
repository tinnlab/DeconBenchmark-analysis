This repository contains code to reproduce the benchmarking results in the article "Fourteen Years of Cellular Deconvolution: Applications, Benchmark, Methodology, and Challenges".

For installing the standalone R package to run deconvolution methods, visit https://github.com/tinnlab/DeconBenchmark.

To get started, follow the steps below:

1. Clone this repository
2. Create a `./data` folder in this repository
3. Download the data from our Zenodo repositories and save them under `./data` folder:
- For Tabula Sapiens: https://zenodo.org/doi/10.5281/zenodo.10687798
- For CELLxGENE: https://zenodo.org/doi/10.5281/zenodo.10688809

**_NOTE 01:_** Users can download the data from CELLxGENE database directly.The following table provides the list of datasets that we used in our analysis. Users can manually download the data by clicking to the link in the `Download URL` column. We also provides the URLs for users to explore the dataset, and the note on choosing the right dataset used in our analysis.

**_NOTE 02:_** Please notice that these data are updated regularly, thus, it may make the repoduced results slightly different to the data deposited on Zenodo.

| Dataset                      | URL                                                                                | Download URL                                                                       | Note                                                  |
|------------------------------|------------------------------------------------------------------------------------|------------------------------------------------------------------------------------|-------------------------------------------------------|
| Tabula Sapiens               | https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5  | https://datasets.cellxgene.cziscience.com/bfd80f12-725c-4482-ad7f-1ed2b4909b0d.rds | Select "Tabula Sapiens - All Cells"                                    |
| Basal Zone of Heart          | https://cellxgene.cziscience.com/collections/43b45a20-a969-49ac-a8e8-8c84b211bd01  | https://datasets.cellxgene.cziscience.com/c3191c1c-e63e-4b9e-b539-25ac669e3c3d.rds |                                                       |
| Fimbria of Uterine Tube      | https://cellxgene.cziscience.com/collections/d36ca85c-3e8b-444c-ba3e-a645040c6185  | https://datasets.cellxgene.cziscience.com/a83b2065-ddb0-41d6-9a46-5e1ed9a9d4ee.rds  | Select "Fallopian tube RNA"                           |
| Primary Visual Cortex        | https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f  | https://datasets.cellxgene.cziscience.com/5b347baa-63d2-4b3e-9d17-5d51c0a3dd23.rds | Select "Dissection: Primary visual cortex(V1)"        |
| Anterior Cingulate Cortex    | https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f  | https://datasets.cellxgene.cziscience.com/5b37a4b8-7bbd-456c-a0ba-40c78411a608.rds | Select "Dissection: Anterior cingulate cortex (ACC)"  |
| Middle Temporal Gyrus        | https://cellxgene.cziscience.com/collections/4dca242c-d302-4dba-a68f-4c61e7bad553  | https://datasets.cellxgene.cziscience.com/fafef79b-a863-4e63-a162-dd933d552ace.rds | Select "Human: Great apes study"                      |
| Primary Auditory Cortex      | https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f  | https://datasets.cellxgene.cziscience.com/4184721b-0a1b-4bb6-b0fb-dc0bf023acb9.rds | Select "Dissection: Primary auditory cortex(A1)"      |
| Primary Somatosensory Cortex | https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f  | https://datasets.cellxgene.cziscience.com/173891cd-830e-4733-861e-5e8ba59973dd.rds | Select "Dissection: Primary somatosensory cortex (S1)" |
| Liver                        | https://cellxgene.cziscience.com/collections/74e10dc4-cbb2-4605-a189-8a1cd8e44d8c  | https://datasets.cellxgene.cziscience.com/ca0b283f-d35d-45d0-91ef-91677bec2663.rds | Select "Lymphoid cells from human liver dataset"      |
| Heart Left Ventricle         | https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb  | https://datasets.cellxgene.cziscience.com/5b6e7a20-b1cb-415e-85cc-5d856ee1591c.rds | Select "DCM/ACM heart cell atlas: All cells" |
| Small Intestine              | https://cellxgene.cziscience.com/collections/7651ac1a-f947-463a-9223-a9e408a41989  | https://datasets.cellxgene.cziscience.com/2c72b4d3-f7b7-4308-8b73-7c5687fd9aa9.rds ||

4. Once the download is finished, each dataset is assigned to a specific ID, which can be seen by running the following command lines in R:

```r
metaFiles <- readRDS("./include/metaFiles.rds")
metaFiles
```

5. Install docker:
   Follow instructions at https://docs.docker.com/engine/install/ to install docker.
   All prebuilt docker images can be found at https://hub.docker.com/u/deconvolution.
   Docker will automatically pull the images when running the scripts below.
   However, CIBERSORT has a specific license requirement and is not available on docker hub.
   Please visit https://github.com/tinnlab/DeconBenchmark-docker to read more about how to build the CIBERSORT docker image
   and any other method.
5Using `R` version >= 4.0.0, run the following commands to install necessary packages used in the analysis:
```r
if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install(c("tidyverse", "babelwhale", "RcppHungarian", "stringr", "processx", "rhdf5", "Seurat", "edgeR"))
```
6. Visit the Matlab license center at http://www.mathworks.com/licensecenter/licenses to obtain a license file.
The hostid of this license is the same as the hostid of the host machine. The user of this license must be root.
This license is required only when running BayesCCE, Deblender, and DecOT.
Save the license file and replace `PATH_TO_MATLAB_LICENSE` in R scripts below with the ABSOLUTE path to the license file.
7. Run `generate-simulated-data-tsp.R` to generate simulated data using single cell data from Tabula Sapiens.
All simulated data will be saved in:
   - `data/sim-tsp` for accuracy and consistency benchmarking
   - `data/scalability` for scalability benchmarking
8. Run `generate-simulated-data-cxg.R` to generate simulated data using single cell data from CELLxGENE.
All simulated data will be saved in:
   - `data/sim-cxg` for accuracy benchmarking
9. Run `accuracy-tsp.R` and `accuracy-cxg.R` to generate the accuracy results.
All results will be saved in the `results/accuracy-tsp` and `results/accuracy-cxg` directories, respectively, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>.rds`.
   Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`. 
   This `rds` file contains a dataframe with method name, dataset name, and 4 accuracy metrics.
10. Run `consistency.R` to generate the consistency results.
   All results will be saved in the `results/consistency` directory, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>_<seed>.rds`.
     Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`.
     This `rds` file contains a dataframe with method name, dataset name, sample, seed, MAE, and correlation.
11. Run `scalability.R` to generate the scalability results.
At the same time, run `memory-monitor.R` to monitor the memory usage of the running docker containers.
All scalability results will be saved in the `results/scalability` directory, including:
   - Individual predicted proportion results for each method and each tissue with format name `<method>_<tissue>_<numberOfBulkSample>.rds`.
     Each `rds` file contains a dataframe of predicted proportions, the ground truth proportions, and the running time.
   - All combined accuracies result with name `allRes.rds`.
     This `rds` file contains a dataframe with method name, dataset name, number of bulk samples, running time, and memory usage.
12. For missing cell types analysis, run `generate-simulated-data-cxg-missingCT.R` to generate simulated data using single cell data from CELLxGENE for the analysis of missing cell types in the reference.
All simulated data will be saved in:
   - `data/sim-cxg-missing25` for accuracy benchmarking when the reference data have 25% of cell types are missing
   - `data/sim-cxg-missing50` for accuracy benchmarking when the reference data have 50% of cell types are missing
13. Run `missing-celltypes-analysis-25.R` and `missing-celltypes-analysis-50.R` to generate the accuracy results in the two scenarios of missing cell types.

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
 [1] edgeR_3.36.0       limma_3.50.3       RcppHungarian_0.2  rhdf5_2.38.1       babelwhale_1.1.0   SeuratObject_4.1.3
 [7] Seurat_4.3.0       lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0        purrr_1.0.1       
[13] readr_2.1.4        tidyr_1.3.0        tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0   

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
