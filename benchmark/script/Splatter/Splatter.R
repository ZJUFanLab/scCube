library(splatter)
library(SingleCellExperiment)
method_list <- c('splat', 'simplesplat', 'kersplat', 'zinb')
###################### DLPFC ######################
slice_list <- c('151507', '151508', '151509', '151510', '151669', '151670', 
                '151671', '151672', '151673', '151674', '151675', '151676')
running_time <- data.frame(type = c('train', 'generate'), 
                           method =  rep(rep(method_list, each=2), 12),
                           time = NA, 
                           slice = rep(slice_list, each=8))


for (i in 1:length(slice_list)) {
  print(paste0('slice: ', slice_list[i]))
  sc_data <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/DLPFC/processed/', slice_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  for (j in 1:length(method_list)) {
    print(paste0('method: ', method_list[j]))
    data <- sc_data[rowSums(sc_data) > 0, ]
    meta <- sc_meta
    sce <- SingleCellExperiment(assays=list(counts = as.matrix(data)), 
                                colData = meta, 
                                rowData = rownames(data))
    set.seed(123)
    

    if(method_list[j] == 'splat'){
      # trian
      t1 <- proc.time()
      params <- splatEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- splatSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'simplesplat'){
      # trian
      t1 <- proc.time()
      params <- simpleEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- simpleSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'kersplat'){
      # trian
      t1 <- proc.time()
      params <- kersplatEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- kersplatSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'zinb'){
      # trian
      t1 <- proc.time()
      params <- zinbEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- zinbSimulate(params)
      t4 <- proc.time()
    }
    
    running_time[running_time$slice == slice_list[i] & 
                   running_time$type == 'train' & 
                   running_time$method == method_list[j], ]$time <- as.numeric(t2 - t1)[3]
    
    running_time[running_time$slice == slice_list[i] & 
                   running_time$type == 'generate' & 
                   running_time$method == method_list[j], ]$time <- as.numeric(t4 - t3)[3]
    
    generate_data <- data.frame(counts(sim))
    rownames(generate_data) <- rownames(sc_data)
    generate_meta <- data.frame(sim@colData)
    generate_meta$Cell_type <- sc_meta$Cell_type
    generate_meta$x <- sc_meta$x
    generate_meta$y <- sc_meta$y
    
    save(generate_data, generate_meta, file = paste0('result/Splatter/', method_list[j], '/', method_list[j], '_DLPFC_', slice_list[i], '_data.Rdata'))
  }
}

write.csv(running_time, file = 'result/Splatter/splatter_DLPFC_self_time.csv')


###################### MERFISH ######################
bregma_list <- c('n0.29', 'n0.24', 'n0.19', 'n0.14', 'n0.09', 'n0.04', '0.01', '0.06', '0.11', '0.16', '0.21', '0.26')
running_time <- data.frame(type = c('train', 'generate'), 
                           method =  rep(rep(method_list, each=2), 12),
                           time = NA, 
                           slice = rep(bregma_list, each=8))


for (i in 1:length(bregma_list)) {
  print(paste0('bregma: ', bregma_list[i]))
  sc_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_data.csv'), row.names = 1)
  sc_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', bregma_list[i], '_meta.csv'), row.names = 1)
  colnames(sc_data) <- rownames(sc_meta)
  
  for (j in 1:length(method_list)) {
    print(paste0('method: ', method_list[j]))
    data <- sc_data[rowSums(sc_data) > 0, ]
    meta <- sc_meta
    if(method_list[j] == 'zinb'){
      data <- round(data)
      data <- data[rowSums(data) > 0, ]
      sce <- SingleCellExperiment(assays=list(counts = as.matrix(data)), 
                                  colData = meta, 
                                  rowData = rownames(data))
    } else{
      sce <- SingleCellExperiment(assays=list(counts = as.matrix(data)), 
                                  colData = meta, 
                                  rowData = rownames(data))
    }
    
    set.seed(123)
    
    if(method_list[j] == 'splat'){
      # trian
      t1 <- proc.time()
      params <- splatEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- splatSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'simplesplat'){
      # trian
      t1 <- proc.time()
      params <- simpleEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- simpleSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'kersplat'){
      # trian
      t1 <- proc.time()
      params <- kersplatEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- kersplatSimulate(params)
      t4 <- proc.time()
    } else if(method_list[j] == 'zinb'){
      # trian
      t1 <- proc.time()
      params <- zinbEstimate(sce)
      t2 <- proc.time()
      
      # generate
      t3 <- proc.time()
      sim <- zinbSimulate(params)
      t4 <- proc.time()
    }
    
    running_time[running_time$slice == bregma_list[i] & 
                   running_time$type == 'train' & 
                   running_time$method == method_list[j], ]$time <- as.numeric(t2 - t1)[3]
    
    running_time[running_time$slice == bregma_list[i] & 
                   running_time$type == 'generate' & 
                   running_time$method == method_list[j], ]$time <- as.numeric(t4 - t3)[3]
    
    generate_data <- data.frame(counts(sim))
    rownames(generate_data) <- rownames(data)
    generate_meta <- data.frame(sim@colData)
    generate_meta$Cell_type <- sc_meta$Cell_type
    generate_meta$x <- sc_meta$x
    generate_meta$y <- sc_meta$y
    
    save(generate_data, generate_meta, file = paste0('result/Splatter/', method_list[j], '/', method_list[j], '_MERFISH_', bregma_list[i], '_data.Rdata'))
  }
}

write.csv(running_time, file = 'result/Splatter/splatter_MERFISH_self_time.csv')

