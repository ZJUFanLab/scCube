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

########################## load evaluate result (DLPFC) ##########################
load("~/workspace/scCube/evaluate/sccube_DLPFC_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/scdesign2_DLPFC_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/splatter_DLPFC_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/symsim_DLPFC_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/srtsim_DLPFC_evaluate_result.Rdata")

# all
all_res <- na.omit(rbind(sccube_res, srtsim_res, symsim_res, scdesign2_res, splatter_res))
all_res <- data.frame(aggregate(all_res$pearson, by=list(all_res$slice, all_res$method), mean))
colnames(all_res) <- c('slice', 'method', 'pcc')
all_res$method <- factor(all_res$method, levels = c('scCube', 'SRTsim', 'scDesign2', 'zinb', 'SymSim', 'splat', 'kersplat', 'simplesplat'))
p <- ggplot(all_res, aes(method, pcc, fill = method)) +
  geom_boxplot() +
  theme_classic()
p <- p + ggbreak::scale_y_break(c(0.015, 0.5),
                                space = 1,
                                scales = 1)
p
ggsave(filename = 'figures/DLPFC_self_benckmark.pdf', p, width = 6, height = 7)


#genes
gene_list <- c('MGP', 'HPCAL1', 'HOPX', 'NEFH', 'PCP4', 'KRT17', 'MOBP')

# ground truth
real_meta <- read.csv(paste0('data/DLPFC/processed/', '151507', '_meta.csv'), row.names = 1)
real_data <- read.csv(paste0('data/DLPFC/processed/', '151507', '_data.csv'), row.names = 1)
rownames(real_meta) <- colnames(real_data)
real_data <- get_normalized_data(data = real_data, meta = real_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = real_data, real_meta = real_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_groundtruth_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# scCube
generate_data <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_', '151507', '_epoch100000_data.csv'), row.names = 1)
generate_meta <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_DLPFC_', '151507', '_epoch100000_meta.csv'), row.names = 1)
rownames(generate_meta) <- colnames(generate_data)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_scCube_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# SRTsim
load(paste0('~/workspace/scCube/result/SRTsim/srtsim_DLPFC_', '151507', '_data.Rdata'))
generate_meta <- data.frame(simSRT1@simcolData)
generate_data <- data.frame(simSRT1@simCounts)
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)

for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_SRTsim_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# scDesign2
load(paste0('~/workspace/scCube/result/scDesign2/scdesign2_DLPFC_', '151507', '_data.Rdata'))
rownames(sc_meta_generate) <- colnames(sim_count_copula_tmp)
# normalize
generate_data <- get_normalized_data(data = sim_count_copula_tmp, meta = sc_meta_generate)

real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
sc_meta_generate$spot <- 'unassigned'
for (j in 1:nrow(sc_meta_generate)) {
  sc_meta_generate[j, ]$spot <- real_meta[real_meta$x == sc_meta_generate[j, ]$x & real_meta$y == sc_meta_generate[j, ]$y, ]$spot
}
rownames(sc_meta_generate) <- colnames(generate_data) <- sc_meta_generate$spot
sc_meta_generate <- sc_meta_generate[real_meta$spot, ]
generate_data <- generate_data[, real_meta$spot]
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_scDesign2_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# zinb
load(paste0('~/workspace/scCube/result/Splatter/', 'zinb', '/', 'zinb', '_DLPFC_', '151507', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_zinb_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# splat
load(paste0('~/workspace/scCube/result/Splatter/', 'splat', '/', 'splat', '_DLPFC_', '151507', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_splat_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# simplesplat
load(paste0('~/workspace/scCube/result/Splatter/', 'simplesplat', '/', 'simplesplat', '_DLPFC_', '151507', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_simplesplat_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

# kersplat
load(paste0('~/workspace/scCube/result/Splatter/', 'kersplat', '/', 'kersplat', '_DLPFC_', '151507', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_kersplat_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}

#symsim
load(paste0('~/workspace/scCube/result/SymSim/symsim_DLPFC_', '151507', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/DLPFC_symsim_', gene_list[i], '.pdf'), p, width = 8, height = 6)
}





########################## load evaluate result (MERFISH) ##########################
load("~/workspace/scCube/evaluate/sccube_MERFISH_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/scdesign2_MERFISH_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/splatter_MERFISH_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/symsim_MERFISH_evaluate_result.Rdata")
load("~/workspace/scCube/evaluate/srtsim_MERFISH_evaluate_result.Rdata")

# all
all_res <- na.omit(rbind(sccube_res, srtsim_res, symsim_res, scdesign2_res, splatter_res))
all_res <- data.frame(aggregate(all_res$pearson, by=list(all_res$slice, all_res$method), mean))
colnames(all_res) <- c('slice', 'method', 'pcc')
all_res$method <- factor(all_res$method, levels = c('scCube', 'SRTsim', 'scDesign2', 'zinb', 'SymSim', 'splat', 'kersplat', 'simplesplat'))
p <- ggplot(all_res, aes(method, pcc, fill = method)) +
  geom_boxplot() +
  theme_classic()
p <- p + scale_y_break(c(0.2, 0.9),
                       space = 1,
                       scales = 1)
p
ggsave(filename = 'figures/MERFISH_self_benckmark.pdf', p, width = 6, height = 7)

# spatial expression patterns of genes (bregma: 0.06 as example)

# OD Mature: Mbp; OD Immature: Pdgfra; Inhibitory: Gad1; Excitatory: Slc17a6; Microglia: Selplg; Astrocyte: Ttyh2; Endothelial: Fn1; Pericytes: Myh11; 
# Ependymal: Nnat
gene_list <- c('Mbp', 'Pdgfra', 'Gad1', 'Slc17a6', 'Selplg', 'Ttyh2', 'Fn1', 'Myh11', 'Nnat', 'Aqp4')


# ground truth
real_meta <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', '0.06', '_meta.csv'), row.names = 1)
real_data <- read.csv(paste0('data/MERFISH/processed/Animal1_Bregma_', '0.06', '_data.csv'), row.names = 1)
rownames(real_meta) <- colnames(real_data)
real_data <- get_normalized_data(data = real_data, meta = real_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = real_data, real_meta = real_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_groundtruth_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# scCube
generate_data <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_MERFISH_', '0.06', '_epoch50000_data.csv'), row.names = 1)
generate_meta <- read.csv(paste0('~/workspace/scCube/result/scCube/sccube_MERFISH_', '0.06', '_epoch50000_meta.csv'), row.names = 1)
rownames(generate_meta) <- colnames(generate_data)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_scCube_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# SRTsim
load(paste0('~/workspace/scCube/result/SRTsim/srtsim_MERFISH_', '0.06', '_data.Rdata'))
generate_meta <- data.frame(simSRT1@simcolData)
generate_data <- data.frame(simSRT1@simCounts)
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)

for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_SRTsim_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# scDesign2
load(paste0('~/workspace/scCube/result/scDesign2/scdesign2_MERFISH_', '0.06', '_data.Rdata'))
rownames(sc_meta_generate) <- colnames(sim_count_copula_tmp)
# normalize
generate_data <- get_normalized_data(data = sim_count_copula_tmp, meta = sc_meta_generate)

real_meta$spot <- paste0('spot_', 1:nrow(real_meta))
sc_meta_generate$spot <- 'unassigned'
for (j in 1:nrow(sc_meta_generate)) {
  sc_meta_generate[j, ]$spot <- real_meta[real_meta$x == sc_meta_generate[j, ]$x & real_meta$y == sc_meta_generate[j, ]$y, ]$spot
}
rownames(sc_meta_generate) <- colnames(generate_data) <- sc_meta_generate$spot
sc_meta_generate <- sc_meta_generate[real_meta$spot, ]
generate_data <- generate_data[, real_meta$spot]
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_scDesign2_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# zinb
load(paste0('~/workspace/scCube/result/Splatter/', 'zinb', '/', 'zinb', '_MERFISH_', '0.06', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_zinb_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# splat
load(paste0('~/workspace/scCube/result/Splatter/', 'splat', '/', 'splat', '_MERFISH_', '0.06', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_splat_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# simplesplat
load(paste0('~/workspace/scCube/result/Splatter/', 'simplesplat', '/', 'simplesplat', '_MERFISH_', '0.06', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_simplesplat_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

# kersplat
load(paste0('~/workspace/scCube/result/Splatter/', 'kersplat', '/', 'kersplat', '_MERFISH_', '0.06', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_kersplat_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}

#symsim
load(paste0('~/workspace/scCube/result/SymSim/symsim_MERFISH_', '0.06', '_data.Rdata'))
rownames(generate_meta) <- colnames(generate_data)
generate_data <- get_normalized_data(data = generate_data, meta = generate_meta)
for (i in 1:length(gene_list)) {
  p <- plot_exp_pattern(real_data = generate_data, real_meta = generate_meta, gene = gene_list[i])
  ggsave(filename = paste0('figures/MERFISH_symsim_', gene_list[i], '.pdf'), p, width = 7.5, height = 6)
}









