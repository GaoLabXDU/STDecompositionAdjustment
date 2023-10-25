# Decom_SPOTlight

if(FALSE){
  "
  The decomposition of SPOTlight methods.
  We downloaded SPOTlight v0.1.7 from https://github.com/MarcElosua/SPOTlight.
  reference:
  [19] Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Res. 2021;49(9):e50. doi:10.1093/nar/gkab043.
  "
}


# SPOTlight ####
Decom_SPOTlight <- function(sc_exp, sc_meta, st_exp, st_xy){
  
  if (!is.data.frame(sc_exp)) stop('ERROR: sc_exp must be a DataFrame!')
  if (!is.data.frame(sc_meta)) stop('ERROR: sc_meta must be a DataFrame!')
  if (!is.data.frame(st_exp)) stop('ERROR: st_exp must be a DataFrame!')
  if (!is.data.frame(st_xy)) stop('ERROR: st_xy must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  
  print('Start decomposition by SPOTlight method.')
  library(dplyr)
  library(SPOTlight)
  library(Seurat)
  library(corrplot)
  
  Celltypes <- unique(sc_meta$celltype)
  
  sc_Seurat <- CreateSeuratObject(counts = sc_exp)
  sc_Seurat@meta.data$celltype <- sc_meta$celltype
  sc_Seurat %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30) %>%
    Seurat::FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.3, verbose = FALSE) -> sc_Seurat
  
  st_Seurat <- CreateSeuratObject(counts = st_exp)
  st_Seurat %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:50) %>%
    Seurat::FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.5, verbose = FALSE) -> st_Seurat
  
  Seurat::Idents(object = sc_Seurat) <- sc_Seurat@meta.data$celltype
  cluster_markers_all <- Seurat::FindAllMarkers(sc_Seurat, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
  
  
  set.seed(123)
  CTprop_tmp <- spotlight_deconvolution(se_sc = sc_Seurat, counts_spatial = st_Seurat@assays$RNA@counts,
                                        clust_vr = 'celltype', cluster_markers = cluster_markers_all, cl_n = 50, hvg = 3000, 
                                        ntop = NULL, transf = 'uv', method = "nsNMF", min_cont = 0.09)
  
  CTprop_tmp1 <- CTprop_tmp[[2]][, 1:length(Celltypes)]
  CTprop_tmp1 <- as.data.frame(CTprop_tmp1)
  rownames(CTprop_tmp1) <- rownames(st_exp)
  SPOTlight_CTprop <- CTprop_tmp1
  return(SPOTlight_CTprop)
}




data_path <- '../Datasets/Simulated data I/'

sc_exp <- read.csv(paste0(data_path, 'simulated st/SC_exp.csv'), header = T, row.names = 1, check.names = F)
sc_meta <- read.csv(paste0(data_path, 'simulated st/SC_meta.csv'), header = T, row.names = 1, check.names = F)
st_exp <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_exp.csv'), header = T, row.names = 1, check.names = F)
st_xy <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_xy.csv'), header = T, row.names = 1, check.names = F)

SPOTlight_CTprop <- Decom_SPOTlight(sc_exp, sc_meta, st_exp, st_xy)


