library(progress)
library(parallel)
library(dplyr)
library(foreach)
########################## function ##########################
cal_pcc <- function(n.core = NULL,
                    real_data,
                    generate_data,
                    slice,
                    method){
  # clear register
  if(is.null(n.core)){
    n.core <- detectCores() - 1
  }
  if (is.null(n.core) || 1 < n.core) {
    # clear register
    doParallel::stopImplicitCluster()
  }
  # register parallel
  if (is.null(n.core)) {
    n.core <- round(parallel::detectCores() / 2)
  }
  n.core <- max(1, n.core)
  if (n.core != 1) {
    doParallel::registerDoParallel(cores = n.core)
  }
  
  res_tmp <- foreach::foreach(n = 1:nrow(real_data), .combine='rbind') %dopar% {
    a <- as.numeric(real_data[n, ])
    b <- as.numeric(generate_data[n, ])
    cor1 <- cor(a, b, method = 'pearson')
    cor1 <- data.frame(pearson = cor1)
    return(cor1)
  }
  
  res_tmp$gene <- rownames(real_data)
  res_tmp$slice <- slice
  res_tmp$method <- method
  
  return(res_tmp)
}


get_normalized_data <- function(data,
                                meta){
  # normalize
  seu <- CreateSeuratObject(data, meta.data = meta)
  seu <- NormalizeData(seu)
  data <- data.frame(seu@assays$RNA@data)
  colnames(data) <- rownames(meta)
  return(data)
}



########################## parameter ##########################
slice_list <- c('151507', '151508', '151509', '151510', '151669', '151670', 
                '151671', '151672', '151673', '151674', '151675', '151676')

bregma_list <- c('n0.29', 'n0.24', 'n0.19', 'n0.14', 'n0.09', 'n0.04', '0.01', '0.06', '0.11', '0.16', '0.21', '0.26')


#############################################################
#                           DLPFC                           #
#############################################################

########################## SRTsim ##########################
srtsim_res <- data.frame()
for (i in 1:length(slice_list)) {
  load(paste0('~/workspace/scCube/result/SRTsim/srtsim_DLPFC_', slice_list[i], '_data.Rdata'))
  generate_meta <- data.frame(simSRT1@simcolData)
  generate_data <- data.frame(simSRT1@simCounts)
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$Cell == generate_meta$Cell)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = slice_list[i],
                       method = 'SRTsim')
  }
  
  srtsim_res <- rbind(srtsim_res, pcc_tmp)
}

save(srtsim_res, file = 'evaluate/srtsim_DLPFC_evaluate_result.Rdata')


########################## SymSim ##########################
symsim_res <- data.frame()
for (i in 1:length(slice_list)) {
  load(paste0('~/workspace/scCube/result/SymSim/symsim_DLPFC_', slice_list[i], '_data.Rdata'))
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$Cell == generate_meta$Cell)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = slice_list[i],
                       method = 'SymSim')
  }
  
  symsim_res <- rbind(symsim_res, pcc_tmp)
}

save(symsim_res, file = 'evaluate/symsim_DLPFC_evaluate_result.Rdata')


########################## Splatter ##########################
splatter_res <- data.frame()
model_list <- c('splat', 'simplesplat', 'kersplat', 'zinb')
for (i in 1:length(slice_list)) {
  real_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  # normalize
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  for (j in 1:length(model_list)) {
    load(paste0('~/workspace/scCube/result/Splatter/', model_list[j], '/', model_list[j], '_DLPFC_', slice_list[i], '_data.Rdata'))
    rownames(generate_meta) <- colnames(generate_data)
    generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
    
    if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == generate_meta$x) & all(real_meta$y == generate_meta$y)){
      pcc_tmp <- cal_pcc(real_data = real_data,
                         generate_data = generate_data,
                         slice = slice_list[i],
                         method = model_list[j])
    }
    
    splatter_res <- rbind(splatter_res, pcc_tmp)
  }
}

save(splatter_res, file = 'evaluate/splatter_DLPFC_evaluate_result.Rdata')


########################## scDesign2 ##########################
scdesign2_res <- data.frame()
for (i in 1:length(slice_list)) {
  load(paste0('~/workspace/scCube/result/scDesign2/scdesign2_DLPFC_', slice_list[i], '_data.Rdata'))
  rownames(sc_meta_generate) <- colnames(sim_count_copula_tmp)
  
  real_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = sim_count_copula_tmp, meta = sc_meta_generate)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
  sc_meta_generate$spot <- 'unassigned'
  for (j in 1:nrow(sc_meta_generate)) {
    sc_meta_generate[j, ]$spot <- real_meta[real_meta$x == sc_meta_generate[j, ]$x & real_meta$y == sc_meta_generate[j, ]$y, ]$spot
  }
  rownames(sc_meta_generate) <- colnames(generate_data) <- sc_meta_generate$spot
  sc_meta_generate <- sc_meta_generate[real_meta$spot, ]
  generate_data <- generate_data[, real_meta$spot]
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == sc_meta_generate$x) & all(real_meta$y == sc_meta_generate$y)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = slice_list[i],
                       method = 'scDesign2')
  }
  
  scdesign2_res <- rbind(scdesign2_res, pcc_tmp)
}

save(scdesign2_res, file = 'evaluate/scdesign2_DLPFC_evaluate_result.Rdata')


########################## scCube ##########################
sccube_res <- data.frame()
for (i in 1:length(slice_list)) {
  generate_data <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_', slice_list[i], '_epoch100000_data.csv'), row.names = 1)
  generate_meta <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_', slice_list[i], '_epoch100000_meta.csv'), row.names = 1)
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
  generate_meta$spot <- 'unassigned'
  
  for (j in 1:nrow(generate_meta)) {
    generate_meta[j, ]$spot <- real_meta[real_meta$x == generate_meta[j, ]$x & real_meta$y == generate_meta[j, ]$y, ]$spot
  }
  rownames(generate_meta) <- colnames(generate_data) <- generate_meta$spot
  generate_meta <- generate_meta[real_meta$spot, ]
  generate_data <- generate_data[, real_meta$spot]
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == generate_meta$x) & all(real_meta$y == generate_meta$y)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = slice_list[i],
                       method ='scCube')
  }
  
  sccube_res <- rbind(sccube_res, pcc_tmp)
}

save(sccube_res, file = 'evaluate/sccube_DLPFC_evaluate_result.Rdata')



#############################################################
#                          MERFISH                         #
#############################################################

########################## SRTsim ##########################
srtsim_res <- data.frame()
for (i in 1:length(bregma_list)) {
  load(paste0('~/workspace/scCube/result/SRTsim/srtsim_MERFISH_', bregma_list[i], '_data.Rdata'))
  generate_meta <- data.frame(simSRT1@simcolData)
  generate_data <- data.frame(simSRT1@simCounts)
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$Cell == generate_meta$Cell)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = bregma_list[i],
                       method = 'SRTsim')
  }
  
  srtsim_res <- rbind(srtsim_res, pcc_tmp)
}

save(srtsim_res, file = 'evaluate/srtsim_MERFISH_evaluate_result.Rdata')


########################## SymSim ##########################
symsim_res <- data.frame()
for (i in 1:length(bregma_list)) {
  load(paste0('~/workspace/scCube/result/SymSim/symsim_MERFISH_', bregma_list[i], '_data.Rdata'))
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$Cell == generate_meta$Cell)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = bregma_list[i],
                       method = 'SymSim')
  }
  
  symsim_res <- rbind(symsim_res, pcc_tmp)
}

save(symsim_res, file = 'evaluate/symsim_MERFISH_evaluate_result.Rdata')


########################## Splatter ##########################
splatter_res <- data.frame()
model_list <- c('splat', 'simplesplat', 'kersplat', 'zinb')
for (i in 1:length(bregma_list)) {
  real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  # normalize
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  for (j in 1:length(model_list)) {
    load(paste0('~/workspace/scCube/result/Splatter/', model_list[j], '/', model_list[j], '_MERFISH_', bregma_list[i], '_data.Rdata'))
    rownames(generate_meta) <- colnames(generate_data)
    generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
    
    if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == generate_meta$x) & all(real_meta$y == generate_meta$y)){
      pcc_tmp <- cal_pcc(real_data = real_data,
                         generate_data = generate_data,
                         slice = bregma_list[i],
                         method = model_list[j])
    }
    
    splatter_res <- rbind(splatter_res, pcc_tmp)
  }
}

save(splatter_res, file = 'evaluate/splatter_MERFISH_evaluate_result.Rdata')


########################## scDesign2 ##########################
scdesign2_res <- data.frame()
for (i in 1:length(bregma_list)) {
  load(paste0('~/workspace/scCube/result/scDesign2/scdesign2_MERFISH_', bregma_list[i], '_data.Rdata'))
  rownames(sc_meta_generate) <- colnames(sim_count_copula_tmp)
  
  real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  # normalize
  generate_data <- get_normalized_data(data = sim_count_copula_tmp, meta = sc_meta_generate)
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
  sc_meta_generate$spot <- 'unassigned'
  for (j in 1:nrow(sc_meta_generate)) {
    sc_meta_generate[j, ]$spot <- real_meta[real_meta$x == sc_meta_generate[j, ]$x & real_meta$y == sc_meta_generate[j, ]$y, ]$spot
  }
  rownames(sc_meta_generate) <- colnames(generate_data) <- sc_meta_generate$spot
  sc_meta_generate <- sc_meta_generate[real_meta$spot, ]
  generate_data <- generate_data[, real_meta$spot]
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == sc_meta_generate$x) & all(real_meta$y == sc_meta_generate$y)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = bregma_list[i],
                       method = 'scDesign2')
  }
  
  scdesign2_res <- rbind(scdesign2_res, pcc_tmp)
}

save(scdesign2_res, file = 'evaluate/scdesign2_MERFISH_evaluate_result.Rdata')


########################## scCube ##########################
sccube_res <- data.frame()
for (i in 1:length(bregma_list)) {
  generate_data <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_MERFISH_', bregma_list[i], '_epoch50000_data.csv'), row.names = 1)
  generate_meta <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_MERFISH_', bregma_list[i], '_epoch50000_meta.csv'), row.names = 1)
  rownames(generate_meta) <- colnames(generate_data)
  
  real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  rownames(real_meta) <- colnames(real_data)
  
  real_data <- get_normalized_data(data = real_data, meta = real_meta)
  
  real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
  generate_meta$spot <- 'unassigned'
  
  for (j in 1:nrow(generate_meta)) {
    generate_meta[j, ]$spot <- real_meta[real_meta$x == generate_meta[j, ]$x & real_meta$y == generate_meta[j, ]$y, ]$spot
  }
  rownames(generate_meta) <- colnames(generate_data) <- generate_meta$spot
  generate_meta <- generate_meta[real_meta$spot, ]
  generate_data <- generate_data[, real_meta$spot]
  
  if(all(rownames(real_data) == rownames(generate_data)) & all(real_meta$x == generate_meta$x) & all(real_meta$y == generate_meta$y)){
    pcc_tmp <- cal_pcc(real_data = real_data,
                       generate_data = generate_data,
                       slice = bregma_list[i],
                       method ='scCube')
  }
  
  sccube_res <- rbind(sccube_res, pcc_tmp)
}

save(sccube_res, file = 'evaluate/sccube_MERFISH_evaluate_result.Rdata')








