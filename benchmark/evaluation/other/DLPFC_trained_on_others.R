########################## function ##########################
plot_exp_pattern <- function(real_data, 
                             real_meta,
                             gene){
  real_meta$gene <- as.numeric(real_data[gene, ])
  real_meta$gene <- (real_meta$gene - min(real_meta$gene)) / (max(real_meta$gene) - min(real_meta$gene))
  p <- ggplot(real_meta, aes(x, y, color = gene)) +
    geom_point(size = 2.5) +
    scale_color_viridis() +
    theme_bw()
  return(p)
}

# train: 151507; test: 151676
# ground truth
real_meta <- read.csv(paste0('data/DLPFC/processed/', '151676', '_meta.csv'), row.names = 1)
real_data <- read.csv(paste0('data/DLPFC/processed/', '151676', '_data.csv'), row.names = 1)
rownames(real_meta) <- colnames(real_data)
real_data <- get_normalized_data(data = real_data, meta = real_meta)
gene_list <- c('MGP', 'HPCAL1', 'HOPX', 'NEFH', 'PCP4', 'KRT17', 'MOBP')
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = real_data, real_meta = real_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_groundtruth_', gene_list[i], '_151676.pdf'), p, width = 8, height = 6)
}


################### SRTsim ###################
library(S4Vectors)
library(SRTsim)
sc_data <- read.csv(paste0('data/DLPFC/processed/', '151507', '_data.csv'), row.names = 1)
sc_meta <- read.csv(paste0('data/DLPFC/processed/', '151507', '_meta.csv'), row.names = 1)
colnames(sc_data) <- rownames(sc_meta)

# fit
example_count <- as.matrix(sc_data)
example_loc <- sc_meta[, c('x', 'y', 'Cell_type')]
colnames(example_loc) <- c("x","y","label")
simSRT <- createSRT(count_in=example_count,loc_in =example_loc)
simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")

test_meta <- read.csv(paste0('data/DLPFC/processed/', '151676', '_meta.csv'), row.names = 1)
example_loc <- test_meta[, c('x', 'y', 'Cell_type')]
colnames(example_loc) <- c("x","y","label")
simSRT1@simcolData <- DataFrame(example_loc)
simSRT1 <- srtsim_count(simSRT1)
save(simSRT1, file = paste0('result/SRTsim/srtsim_DLPFC_151676_data_trained_on_other.Rdata'))

generate_meta <- data.frame(simSRT1@simcolData)
generate_data <- data.frame(simSRT1@simCounts)
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)

pcc_srtsim <- cal_pcc(real_data = real_data,
                      generate_data = generate_data,
                      slice = '151676',
                      method = 'SRTsim')

for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_SRTsim_', gene_list[i], '_151676.pdf'), p, width = 8, height = 6)
}


################### sccube ###################
generate_data <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_151676_epoch100000_data_trained_on_other.csv'), row.names = 1)
generate_meta <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_151676_epoch100000_meta_trained_on_other.csv'), row.names = 1)
rownames(generate_meta) <- colnames(generate_data)
real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
generate_meta$spot <- 'unassigned'

for (j in 1:nrow(generate_meta)) {
  generate_meta[j, ]$spot <- real_meta[real_meta$x == generate_meta[j, ]$x & real_meta$y == generate_meta[j, ]$y, ]$spot
}
rownames(generate_meta) <- colnames(generate_data) <- generate_meta$spot
generate_meta <- generate_meta[real_meta$spot, ]
generate_data <- generate_data[, real_meta$spot]

pcc_sccube <- cal_pcc(real_data = real_data,
                      generate_data = generate_data,
                      slice = '151676',
                      method = 'scCube')

for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_scCube_', gene_list[i], '_151676.pdf'), p, width = 8, height = 6)
}
