library(SymSim)
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
  # fit
  t1 <- proc.time()
  best_params <- BestMatchParams('UMI', sc_data, paste0('slice:', slice_list[i]), n_optimal = 1)
  true_counts_res <- SimulateTrueCounts(ncells_total = ncol(sc_data),
                                        ngenes = nrow(sc_data),
                                        evf_type = 'one.population',
                                        Sigma = best_params$Sigma,
                                        gene_effects_sd = best_params$gene_effects_sd,
                                        gene_effect_prob = best_params$gene_effect_prob,
                                        scale_s = best_params$scale_s,
                                        prop_hge = best_params$prop_hge,
                                        mean_hge = best_params$mean_hge,
                                        randseed = 12345)
  t2 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, nrow(sc_data), replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts = true_counts_res[[1]],
                                         meta_cell = true_counts_res[[3]],
                                         protocol = best_params$protocol,
                                         alpha_mean = best_params$alpha_mean,
                                         alpha_sd = best_params$alpha_sd,
                                         depth_mean = best_params$depth_mean,
                                         depth_sd = best_params$depth_sd,
                                         nPCR1 = best_params$nPCR1,
                                         gene_len = gene_len)
  t4 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  
  generate_data <- data.frame(observed_counts$counts)
  rownames(generate_data) <- rownames(sc_data)
  colnames(generate_data) <- colnames(sc_data)
  generate_meta <- sc_meta
  rownames(generate_meta) <- colnames(generate_data)
  save(generate_data, generate_meta, file = paste0('result/SymSim/symsim_DLPFC_', slice_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/SymSim/symsim_DLPFC_self_time.csv')


###################### MERFISH ######################
bregma_list <- c('n0.29', 'n0.24', 'n0.19', 'n0.14', 'n0.09', 'n0.04', '0.01', '0.06', '0.11', '0.16', '0.21', '0.26')
running_time <- data.frame(type = c('train', 'generate'), 
                           time = NA, 
                           slice = rep(bregma_list, each=2))


for (i in 1:length(bregma_list)) {
  sc_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  # fit
  t1 <- proc.time()
  best_params <- BestMatchParams('UMI', sc_data, paste0('bregma:', bregma_list[i]), n_optimal = 1)
  true_counts_res <- SimulateTrueCounts(ncells_total = ncol(sc_data),
                                        ngenes = nrow(sc_data),
                                        evf_type = 'one.population',
                                        Sigma = best_params$Sigma,
                                        gene_effects_sd = best_params$gene_effects_sd,
                                        gene_effect_prob = best_params$gene_effect_prob,
                                        scale_s = best_params$scale_s,
                                        prop_hge = best_params$prop_hge,
                                        mean_hge = best_params$mean_hge,
                                        randseed = 12345)
  t2 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, nrow(sc_data), replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts = true_counts_res[[1]],
                                         meta_cell = true_counts_res[[3]],
                                         protocol = best_params$protocol,
                                         alpha_mean = best_params$alpha_mean,
                                         alpha_sd = best_params$alpha_sd,
                                         depth_mean = best_params$depth_mean,
                                         depth_sd = best_params$depth_sd,
                                         nPCR1 = best_params$nPCR1,
                                         gene_len = gene_len)
  t4 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  
  generate_data <- data.frame(observed_counts$counts)
  rownames(generate_data) <- rownames(sc_data)
  colnames(generate_data) <- colnames(sc_data)
  generate_meta <- sc_meta
  rownames(generate_meta) <- colnames(generate_data)
  save(generate_data, generate_meta, file = paste0('result/SymSim/symsim_MERFISH_', bregma_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/SymSim/symsim_MERFISH_self_time.csv')


