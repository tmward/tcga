# tcga

This repository holds software used to perform analyses and generate
figures/tables for the 2021 American Society of Colon and Rectal
Surgeons (ASCRS) Annual Scientific Meeting presentation, number SP45,
titled: "Tumor Genotypes Account for Survival Differences in Right and
Left-sided Colon Cancers". This presentation will be presented in the
Plenary Abstract Session III: Colorectal Cancer and Other Neoplasia.
Additionally, it was submitted for publication consideration to Diseases
of the Colon and Rectum.

# Overview of repository contents

## `data` folder

Only contains `gene.csv`, a CSV that contains data from The Cancer
Genome Atlas Program (TCGA) for colorectal cancers (CRC).

## `output` folder

Empty directory that will contain generated R objects, tables, and
figures.

## `src` folder

All scripts expect that they are executed with a working directory of
`src` (e.g., should type `./script.R` to execute, *not* `src/script.R`).

All tables produced are HTML tables. To use in Microsoft Word, it is
best to open them in a browser then copy and paste them into Word.

All plots are in png format. Graphic size and format is easy to tweak by
changing the script.

All file outputs are put into the `../output/` directory with self-explanatory
names.

### `tidy_tcga.R`

Cleans and tidies the CSV holding information on all CRC obtained from
TCGA. Will save a serialized R object (`.rds`) for a tibble holding
tidied information with patient characteristics and genetics and
a separate tibble holding information on MSI. It only keeps patient data
when the patient has information on a primary CRC tumor and completed
genetic mutation data.

### `make_inputs.R`

Produces two R objects from the previously tidied patient and MSI
data produced by `tidy_tcga.R`. Notably, it only keeps data for patients
with a right or left colon cancer and completed information for use in
survival analyses (so more than 0 months of follow-up) and completed T
and N staging. The first, `analyses_data_tidy.rds`, holds the subset of
data for R and L cancers with imputed missing LVI, R0 margin status,
and metastases status. This is tidy (one gene mutation per line). It
also produces a second R object, `model_inputs.rds`, that holds the data
widened, for use in modeling.

### `chisq.R`

Performs the analyses and generates summary tables for them for
chi-squared tests of independence. Requires `tidy_tcga.R` then
`make_inputs.R` to be run first. It analyses:

1. Tumor sidedness versus gene mutational status
2. Nodal positivity (N0 vs N1+) versus gene mutational status

### `survival_analyses.R`

Requires `tidy_tcga.R` then `make_inputs.R` to be run first. This script
will take `model_inputs.rds` and output:

1. Table with lasso regression for 5-yr survival results

### `causal_analyses.R`

Requires `tidy_tcga.R` then `make_inputs.R` to be run first. This script
will take `model_inputs.rds` and output:

1. Table with total effect of R vs L on 5-yr survival Cox PH results
2. Schoenfeld residual plots for the above analysis
3. Table with total effect of genetics on 5-yr survival Cox PH results
4. Schoenfeld residual plots for the above analysis
5. Plot of KM curve of R vs L and 5-yr survival
6. Table with survival statistics each year of follow-up

### `table1.R`

Requires the outputs from `tidy_tcga.R` and `make_inputs.R`. It takes
`model_inputs.rds` and `tidied.rds` and outputs:

1. Table 1 in the paper

# sessionInfo

All analyses were run in the following environment:

```
R version 4.0.3 (2020-10-10)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Fedora 33 (Workstation Edition)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libflexiblas.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] survminer_0.4.8 ggpubr_0.4.0    gtsummary_1.3.5 survival_3.2-7 
 [5] glmnet_4.0-2    Matrix_1.2-18   forcats_0.5.0   stringr_1.4.0  
 [9] dplyr_1.0.2     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2    
[13] tibble_3.0.4    ggplot2_3.3.2   tidyverse_1.3.0 nvimcom_0.9-102

loaded via a namespace (and not attached):
 [1] fs_1.5.0                 usethis_1.6.3           
 [3] lubridate_1.7.9          httr_1.4.2              
 [5] tools_4.0.3              backports_1.1.10        
 [7] utf8_1.1.4               R6_2.4.1                
 [9] adaptMCMC_1.3            DBI_1.1.0               
[11] colorspace_1.4-1         withr_2.3.0             
[13] tidyselect_1.1.0         gridExtra_2.3           
[15] curl_4.3                 compiler_4.0.3          
[17] cli_2.1.0                rvest_0.3.6             
[19] gt_0.2.2                 xml2_1.3.2              
[21] labeling_0.4.2           sass_0.2.0              
[23] scales_1.1.1             checkmate_2.0.0         
[25] survMisc_0.5.5           commonmark_1.7          
[27] digest_0.6.26            foreign_0.8-80          
[29] rio_0.5.16               pkgconfig_2.0.3         
[31] htmltools_0.5.0          highr_0.8               
[33] dbplyr_1.4.4             rlang_0.4.8             
[35] readxl_1.3.1             selectiveInference_1.2.5
[37] rstudioapi_0.11          shape_1.4.5             
[39] generics_0.0.2           farver_2.0.3            
[41] zoo_1.8-8                jsonlite_1.7.1          
[43] zip_2.1.1                car_3.0-10              
[45] magrittr_1.5             Rcpp_1.0.5              
[47] munsell_0.5.0            fansi_0.4.1             
[49] clipr_0.7.1              abind_1.4-5             
[51] lifecycle_0.2.0          stringi_1.5.3           
[53] carData_3.0-4            MASS_7.3-53             
[55] grid_4.0.3               blob_1.2.1              
[57] parallel_4.0.3           crayon_1.3.4            
[59] lattice_0.20-41          haven_2.3.1             
[61] splines_4.0.3            hms_0.5.3               
[63] knitr_1.30               ps_1.4.0                
[65] pillar_1.4.6             tcltk_4.0.3             
[67] ggsignif_0.6.0           codetools_0.2-16        
[69] reprex_0.3.0             glue_1.4.2              
[71] data.table_1.13.2        modelr_0.1.8            
[73] vctrs_0.3.4              foreach_1.5.1           
[75] cellranger_1.1.0         gtable_0.3.0            
[77] km.ci_0.5-2              assertthat_0.2.1        
[79] xfun_0.18                openxlsx_4.2.3          
[81] xtable_1.8-4             broom_0.7.2             
[83] Rmpfr_0.8-1              coda_0.19-4             
[85] rstatix_0.6.0            intervals_0.15.2        
[87] iterators_1.0.13         KMsurv_0.1-5            
[89] gmp_0.6-1                ellipsis_0.3.1          
```
