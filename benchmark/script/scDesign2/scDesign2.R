library(scDesign2) # copula
###################### DLPFC ######################
slice_list <- c('151507', '151508', '151509', '151510', '151669', '151670', 
                '151671', '151672', '151673', '151674', '151675', '151676')
running_time <- data.frame(type = c('train', 'generate'), 
                           time = NA, 
                           slice = rep(slice_list, each=2))

for (i in 1:length(slice_list)) {
  sc_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  # fit model and simulate data -----------------------------------------------------------
  set.seed(1)
  sc_data2 <- sc_data
  colnames(sc_data2) <- rep('cell', ncol(sc_data2))
  t1 <- proc.time()
  copula_result <- fit_model_scDesign2(as.matrix(sc_data2), 'cell', sim_method = 'copula')
  t2 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  sim_count_copula_tmp <- simulate_count_scDesign2(copula_result, 
                                                   n_cell_new = copula_result$cell$n_cell, 
                                                   sim_method = 'copula')
  t4 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  sc_meta_generate <- sc_meta
  colnames(sim_count_copula_tmp) <- sc_meta_generate$Cell
  rownames(sim_count_copula_tmp) <- rownames(sc_data2)
  
  save(sim_count_copula_tmp, sc_meta_generate, file = paste0('result/scDesign2/scdesign2_DLPFC_', slice_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/scDesign2/scdesign2_DLPFC_self_time.csv')


###################### MERFISH ######################
bregma_list <- c('n0.29', 'n0.24', 'n0.19', 'n0.14', 'n0.09', 'n0.04', '0.01', '0.06', '0.11', '0.16', '0.21', '0.26')
running_time <- data.frame(type = c('train', 'generate'), 
                           time = NA, 
                           slice = rep(bregma_list, each=2))

for (i in 1:length(bregma_list)) {
  sc_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  # fit model and simulate data -----------------------------------------------------------
  set.seed(1)
  sc_data2 <- sc_data
  colnames(sc_data2) <- rep('cell', ncol(sc_data2))
  t1 <- proc.time()
  copula_result <- fit_model_scDesign2(as.matrix(sc_data2), 'cell', sim_method = 'copula')
  t2 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  sim_count_copula_tmp <- simulate_count_scDesign2(copula_result, 
                                                   n_cell_new = copula_result$cell$n_cell, 
                                                   sim_method = 'copula')
  t4 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  sc_meta_generate <- sc_meta
  colnames(sim_count_copula_tmp) <- sc_meta_generate$Cell
  
  save(sim_count_copula_tmp, sc_meta_generate, file = paste0('result/scDesign2/scdesign2_MERFISH_', bregma_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/scDesign2/scdesign2_MERFISH_self_time.csv')



