library(S4Vectors)
library(SRTsim)
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
  example_count <- as.matrix(sc_data)
  example_loc <- sc_meta[, c('x', 'y', 'Cell_type')]
  colnames(example_loc) <- c("x","y","label")
  t1 <- proc.time()
  simSRT <- createSRT(count_in=example_count,loc_in =example_loc)
  set.seed(123)
  simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
  t2 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  simSRT1 <- srtsim_count(simSRT1)
  t4 <- proc.time()
  running_time[running_time$slice == slice_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  save(simSRT1, file = paste0('result/SRTsim/srtsim_DLPFC_', slice_list[i], '_data.Rdata'))
}

write.csv(running_time, file = 'result/SRTsim/srtsim_DLPFC_self_time.csv')

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
  example_count <- as.matrix(sc_data)
  example_loc <- sc_meta[, c('x', 'y', 'Cell_type')]
  colnames(example_loc) <- c("x","y","label")
  t1 <- proc.time()
  simSRT <- createSRT(count_in=example_count,loc_in =example_loc)
  set.seed(123)
  simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
  t2 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'train', ]$time <- as.numeric(t2 - t1)[3]
  
  # generate
  t3 <- proc.time()
  simSRT1 <- srtsim_count(simSRT1)
  t4 <- proc.time()
  running_time[running_time$slice == bregma_list[i] & running_time$type == 'generate', ]$time <- as.numeric(t4 - t3)[3]
  
  save(simSRT1, file = paste0('result/SRTsim/srtsim_MERFISH_', bregma_list[i], '_data.Rdata'))
}
write.csv(running_time, file = 'result/SRTsim/srtsim_MERFISH_self_time.csv')
