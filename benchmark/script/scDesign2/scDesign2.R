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
  colnames(sc_data2) <- sc_meta$Cell_type
  t1 <- proc.time()
  copula_result <- fit_model_scDesign2(as.matrix(sc_data2), unique(sc_meta$Cell_type), sim_method = 'copula')
  t2 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  sc_meta$Cell_type <- factor(sc_meta$Cell_type, levels = unique(sc_meta$Cell_type))
  t3 <- proc.time()
  sim_count_copula_tmp <- simulate_count_scDesign2(copula_result, 
                                                   n_cell_new = nrow(sc_meta), 
                                                   cell_type_prop = as.numeric(table(sc_meta$Cell_type))/nrow(sc_meta), 
                                                   sim_method = 'copula')
  t4 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  # sc_meta_generate <- data.frame(Cell = paste0('C_', 1:ncol(sim_count_copula_tmp)),
  #                         Cell_type = colnames(sim_count_copula_tmp))
  # rownames(sim_count_copula_tmp) <- rownames(sc_data)
  # colnames(sim_count_copula_tmp) <- rownames(sc_meta_generate) <- sc_meta_generate$Cell
  # 
  # sc_meta_generate$x <- 0
  # sc_meta_generate$y <- 0
  # for (j in unique(sc_meta$Cell_type)) {
  #   sc_meta_generate[sc_meta_generate$Cell_type == j, ]$x <- sc_meta[sc_meta$Cell_type == j, ]$x
  #   sc_meta_generate[sc_meta_generate$Cell_type == j, ]$y <- sc_meta[sc_meta$Cell_type == j, ]$y
  # }
  # 
  # save(sim_count_copula_tmp, sc_meta_generate, file = paste0('result/scDesign2/scdesign2_DLPFC_', slice_list[i], '_data.Rdata'))
  
  save(sim_count_copula_tmp, file = paste0('result/scDesign2/scdesign2_DLPFC_', slice_list[i], '_data.Rdata'))
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
  colnames(sc_data2) <- sc_meta$Cell_type
  t1 <- proc.time()
  copula_result <- fit_model_scDesign2(as.matrix(sc_data2), unique(sc_meta$Cell_type), sim_method = 'copula')
  t2 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  sc_meta$Cell_type <- factor(sc_meta$Cell_type, levels = unique(sc_meta$Cell_type))
  t3 <- proc.time()
  sim_count_copula_tmp <- simulate_count_scDesign2(copula_result, 
                                                   n_cell_new = nrow(sc_meta), 
                                                   cell_type_prop = as.numeric(table(sc_meta$Cell_type))/nrow(sc_meta), 
                                                   sim_method = 'copula')
  t4 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  # sc_meta_generate <- data.frame(Cell = paste0('C_', 1:ncol(sim_count_copula_tmp)),
  #                                Cell_type = colnames(sim_count_copula_tmp))
  # rownames(sim_count_copula_tmp) <- rownames(sc_data)
  # colnames(sim_count_copula_tmp) <- rownames(sc_meta_generate) <- sc_meta_generate$Cell
  # 
  # sc_meta_generate$x <- 0
  # sc_meta_generate$y <- 0
  # for (j in unique(sc_meta$Cell_type)) {
  #   sc_meta_generate[sc_meta_generate$Cell_type == j, ]$x <- sc_meta[sc_meta$Cell_type == j, ]$x
  #   sc_meta_generate[sc_meta_generate$Cell_type == j, ]$y <- sc_meta[sc_meta$Cell_type == j, ]$y
  # }
  # 
  # save(sim_count_copula_tmp, sc_meta_generate, file = paste0('result/scDesign2/scdesign2_MERFISH_', bregma_list[i], '_data.Rdata'))
  
  save(sim_count_copula_tmp, file = paste0('result/scDesign2/scdesign2_MERFISH_', bregma_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/scDesign2/scdesign2_MERFISH_self_time.csv')
