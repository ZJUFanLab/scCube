library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
###################### DLPFC ######################
# slice_list <- c('151507', '151508', '151509', '151510', '151669', '151670', 
#                 '151671', '151672', '151673', '151674', '151675', '151676')
slice_list <- c('151507')
running_time <- data.frame(time = NA, slice = slice_list)
# example_sce <- readRDS((url("https://figshare.com/ndownloader/files/40582019")))

for (i in 1:length(slice_list)) {
  sc_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  sce <- SingleCellExperiment(assays=list(counts=sc_data), colData = sc_meta)
  # test
  # sce <- sce[c('GFAP', 'HPCAL1', 'HOPX', 'NEFH', 'PCP4', 'KRT17', 'MOBP'), ]
  set.seed(123)
  t1 <- proc.time()
  example_simu <- scdesign3(
    sce = sce,
    assay_use = "counts",
    celltype = "Cell_type",
    pseudotime = NULL,
    spatial = c("x", "y"),
    other_covariates = NULL,
    mu_formula = "s(x, y, bs = 'gp', k= 400)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 60,
    usebam = FALSE,
    corr_formula = "1",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE
  )
  t2 <- proc.time()
  running_time[running_time$slice == slice_list[i], ]$time <- as.numeric(t2 - t1)[3]
  
  save(example_simu, file = paste0('result/scDesign3/scdesign3_DLPFC_', slice_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/scDesign3/scdesign3_DLPFC_self_time.csv')

###################### MERFISH ######################
bregma_list <- c('n0.29', 'n0.24', 'n0.19', 'n0.14', 'n0.09', 'n0.04', '0.01', '0.06', '0.11', '0.16', '0.21', '0.26')
running_time <- data.frame(time = NA, slice = bregma_list)

for (i in 1:length(bregma_list)) {
  sc_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  sce <- SingleCellExperiment(assays=list(counts=sc_data), colData = sc_meta)
  # test
  set.seed(123)
  t1 <- proc.time()
  example_simu <- scdesign3(
    sce = sce,
    assay_use = "counts",
    celltype = "Cell_type",
    pseudotime = NULL,
    spatial = c("x", "y"),
    other_covariates = NULL,
    mu_formula = "s(x, y, bs = 'gp', k= 400)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 48,
    usebam = FALSE,
    corr_formula = "1",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE
  )
  t2 <- proc.time()
  running_time[running_time$slice == bregma_list[i], ]$time <- as.numeric(t2 - t1)[3]
  
  save(example_simu, file = paste0('result/scDesign3/scdesign3_MERFISH_', bregma_list[i], '_data.Rdata'))
}
write.csv(running_time, file = 'result/scDesign3/scdesign3_MERFISH_self_time.csv')
