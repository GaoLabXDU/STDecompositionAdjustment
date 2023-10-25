# Decom_RCTD

if(FALSE){
  "
  The decomposition of RCTD methods.
  We downloaded RCTD from https://github.com/dmcable/spacexr.
  reference:
  [24] Cable DM, Murray E, Zou LS, Goeva A, Macosko EZ, Chen F, et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat Biotechnol. 2022;40(4):517–526. doi:10.1038/s41587-021-00830-w.
  "
}


# RCTD ####
Decom_RCTD <- function(sc_exp, sc_meta, st_exp, st_xy){
  
  if (!is.data.frame(sc_exp)) stop('ERROR: sc_exp must be a DataFrame!')
  if (!is.data.frame(sc_meta)) stop('ERROR: sc_meta must be a DataFrame!')
  if (!is.data.frame(st_exp)) stop('ERROR: st_exp must be a DataFrame!')
  if (!is.data.frame(st_xy)) stop('ERROR: st_xy must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  
  print('Start decomposition by RCTD method.')
  library(spacexr)
  
  if (!('cellname' %in% colnames(sc_meta))) 
    sc_meta$cellname <- rownames(sc_meta)
  if (!('nUMI' %in% colnames(sc_meta))) {
    sc_nUMI <- colSums(sc_exp)
    sc_nUMI <- as.numeric(sc_nUMI)
    names(sc_nUMI) <- sc_meta$cellname
  } else {
    sc_nUMI <- as.factor(sc_meta$nUMI)
    sc_nUMI <- as.numeric(sc_nUMI)
    names(sc_nUMI) <- sc_meta$cellname
  }
  
  # Create RCTD object
  st_RCTD <- SpatialRNA(coords = st_xy, counts = st_exp, nUMI = colSums(st_exp), require_int = FALSE)
  
  cell_types <- as.factor(sc_meta$celltype)
  names(cell_types) <- sc_meta$cellname
  sc_RCTD <- Reference(counts = sc_exp, cell_types, nUMI = sc_nUMI)
  
  myRCTD <- create.RCTD(st_RCTD, sc_RCTD, max_cores = 1, test_mode = FALSE)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  results <- myRCTD@results
  
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')
  CTprop_tmp1 <- norm_weights@x
  CTprop_tmp2 <- matrix(CTprop_tmp1, nrow = dim(st_exp)[2], ncol = length(unique(cell_types)), 
                        dimnames = list(colnames(st_exp), norm_weights@Dimnames[[2]]))
  RCTD_CTprop <- CTprop_tmp2
  return(RCTD_CTprop)
}





data_path <- '../Datasets/Simulated data I/'

sc_exp <- read.csv(paste0(data_path, 'simulated st/SC_exp.csv'), header = T, row.names = 1, check.names = F)
sc_meta <- read.csv(paste0(data_path, 'simulated st/SC_meta.csv'), header = T, row.names = 1, check.names = F)
st_exp <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_exp.csv'), header = T, row.names = 1, check.names = F)
st_xy <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_xy.csv'), header = T, row.names = 1, check.names = F)

RCTD_CTprop <- Decom_RCTD(sc_exp, sc_meta, st_exp, st_xy)


