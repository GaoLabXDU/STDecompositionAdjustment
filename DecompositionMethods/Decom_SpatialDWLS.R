# Decom_SpatialDWLS

if(FALSE){
  "
  The decomposition of SpatialDWLS methods.
  We downloaded SpatialDWLS in https://github.com/RubD/Giotto and followed the guidelines on https://github.com/rdong08/spatialDWLS_dataset.
  reference:
  [20] Dong R, Yuan G-C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 2021;22(1):145. doi:10.1186/s13059-021-02362-7.
  "
}

# SpatialDWLS ####
Decom_SpatialDWLS <- function(sc_exp, sc_meta, st_exp, st_xy){
  
  if (!is.data.frame(sc_exp)) stop('ERROR: sc_exp must be a DataFrame!')
  if (!is.data.frame(sc_meta)) stop('ERROR: sc_meta must be a DataFrame!')
  if (!is.data.frame(st_exp)) stop('ERROR: st_exp must be a DataFrame!')
  if (!is.data.frame(st_xy)) stop('ERROR: st_xy must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  
  print('Start decomposition by SpatialDWLS method.')
  library(Giotto)
  library(Matrix)
  
  instrs = createGiottoInstructions(python_path = 'D:/anaconda/envs/Giotto-env/python.exe')
  
  # Create Giotto object
  sc_Gio <- createGiottoObject(raw_exprs = sc_exp, instructions = instrs)
  sc_Gio <- normalizeGiotto(gobject = sc_Gio)
  sc_Gio@cell_metadata$leiden_clus <- as.character(sc_meta$celltype)
  gini_markers_subclusters <- findMarkers_one_vs_all(gobject = sc_Gio, method = 'gini', logFC = 0.5, 
                                                     expression_values = 'normalized', cluster_column = 'leiden_clus')
  topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']
  sc_Gio_norm_exp <- 2^(sc_Gio@norm_expr)-1
  
  ExprSubset <- sc_Gio_norm_exp[as.character(topgenes_gini$genes),]
  Sig <- NULL
  for (i in as.character(unique(sc_meta$celltype))){
    Sig <- cbind(Sig, (apply(ExprSubset,1,function(y) mean(y[which(sc_meta$celltype==i)]))))
  }
  colnames(Sig) <- as.character(unique(sc_meta$celltype))
  
  st_Gio <- createGiottoObject(raw_exprs = st_exp, instructions = instrs, spatial_locs = st_xy)
  st_Gio <- normalizeGiotto(gobject = st_Gio)
  st_Gio <- calculateHVG(gobject = st_Gio)
  gene_metadata = fDataDT(st_Gio)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  st_Gio <- runPCA(gobject = st_Gio, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_Gio, genes_to_use = featgenes, scale_unit = F)
  st_Gio <- createNearestNetwork(gobject = st_Gio, dimensions_to_use = 1:10, k = 4)
  st_Gio <- doLeidenCluster(gobject = st_Gio, resolution = 0.4)
  st_Gio <- runDWLSDeconv(gobject = st_Gio, sign_matrix = Sig)
  
  CTprop_tmp <- as.data.frame(st_Gio@spatial_enrichment$DWLS)
  rownames(CTprop_tmp) <- CTprop_tmp[,'cell_ID']
  SpatialDWLS_CTprop <- CTprop_tmp[, -1]
  return(SpatialDWLS_CTprop)
}




data_path <- '../Datasets/Simulated data I/'

sc_exp <- read.csv(paste0(data_path, 'simulated st/SC_exp.csv'), header = T, row.names = 1, check.names = F)
sc_meta <- read.csv(paste0(data_path, 'simulated st/SC_meta.csv'), header = T, row.names = 1, check.names = F)
st_exp <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_exp.csv'), header = T, row.names = 1, check.names = F)
st_xy <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_xy.csv'), header = T, row.names = 1, check.names = F)

SpatialDWLS_CTprop <- Decom_SpatialDWLS(sc_exp, sc_meta, st_exp, st_xy)





