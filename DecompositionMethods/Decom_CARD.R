# Decom_CARD

if(FALSE){
  "
  The decomposition of CARD methods.
  We downloaded CARD from https://github.com/YingMa0107/CARD and followed the instructions in https://yingma0107.github.io/CARD/documentation/04_CARD_Example.html.
  reference:
  [21] Ma Y, Zhou X. Spatially informed cell-type deconvolution for spatial transcriptomics. Nat Biotechnol. 2022;40(9):1349–1359. doi:10.1038/s41587-022-01273-7.
  "
}


# CARD ####
Decom_CARD <- function(sc_exp, sc_meta, st_exp, st_xy){

  if (!is.data.frame(sc_exp)) stop('ERROR: sc_exp must be a DataFrame!')
  if (!is.data.frame(sc_meta)) stop('ERROR: sc_meta must be a DataFrame!')
  if (!is.data.frame(st_exp)) stop('ERROR: st_exp must be a DataFrame!')
  if (!is.data.frame(st_xy)) stop('ERROR: st_xy must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  if (!('X' %in% colnames(st_meta))) stop('ERROR: the colnames of st_meta must have X !')
  if (!('Y' %in% colnames(st_meta))) stop('ERROR: the colnames of st_meta must have Y !')
  
  print('Start decomposition by CARD method.')
  
  library(CARD)
  
  sc_meta$sampleInfo <- rep('sample1', times = dim(sc_meta)[1])
  st_xy <- st_xy %>% rename(x = X, y = X)
  
  # Creat CARD object
  CARD_obj <- createCARDObject(sc_count = as.matrix(sc_exp), sc_meta = sc_meta,
                               spatial_count = as.matrix(st_exp), spatial_location = st_xy,
                               ct.varname = "celltype", ct.select = unique(sc_meta$celltype), 
                               sample.varname = "sampleInfo", minCountGene = 100, minCountSpot = 5) 
  
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  CTprop_tmp <- as.data.frame(CARD_obj@Proportion_CARD[,])
  CARD_CTprop <- CTprop_tmp
  return(CARD_CTprop)
}




data_path <- '../Datasets/Simulated data I/'

sc_exp <- read.csv(paste0(data_path, 'simulated st/SC_exp.csv'), header = T, row.names = 1, check.names = F)
sc_meta <- read.csv(paste0(data_path, 'simulated st/SC_meta.csv'), header = T, row.names = 1, check.names = F)
st_exp <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_exp.csv'), header = T, row.names = 1, check.names = F)
st_xy <- read.csv(paste0(data_path, 'simulated st/ST_Simulated_xy.csv'), header = T, row.names = 1, check.names = F)

CARD_CTprop <- Decom_CARD(sc_exp, sc_meta, st_exp, st_xy)


