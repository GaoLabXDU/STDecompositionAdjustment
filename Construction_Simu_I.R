
#' The construction of Simulated data I

if(FALSE){
  "
  Mouse somatosensory cortex tissue.
  scRNA-seq: 19972 genes in 1691 cells, 6 cell types.
  single cell ST: 10000 genes in 523 cells.
  simulated ST: 10000 genes in 85 spots.
  reference: Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq
  "
}

library(dplyr)
library(Seurat)
library(stringr)
library(Giotto)

#' Simulated data I_Construction Function

# Preprocessing to select cells in ST data and select genes in scRNA-seq data ####
# @st_exp_raw: DataFrame of Raw ST data.
# @st_meta_raw: Raw ST metadata.
# @sc_exp_raw: Raw scRNA-seq data.
# @sc_meta_raw: Raw scRNA-seq metadata.
# @Return: A list of Data_Preprocessed, st_exp, st_meta, sc_exp, sc_meta.

Data_Preprocessing <- function(st_exp_raw, st_meta_raw, sc_exp_raw, sc_meta_raw){
  
  if (!is.data.frame(st_exp_raw)) stop('ERROR: st_exp_raw must be a DataFrame!')
  if (!is.data.frame(st_meta_raw)) stop('ERROR: st_meta_raw must be a DataFrame!')
  if (!is.data.frame(sc_exp_raw)) stop('ERROR: sc_exp_raw must be a DataFrame!')
  if (!is.data.frame(sc_meta_raw)) stop('ERROR: sc_meta_raw must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta_raw)))  
    stop('ERROR: the colnames of sc_meta_raw must have celltype to represent cell type information!')
  
  st_meta_raw <- as.data.frame(st_meta_raw)
  print('Start data preprocessing in Simulated data I.')
  
  # Select cells from 'Cortex' region in st data
  st_cortex_meta <- st_meta_raw[which(st_meta_raw$Region == "Cortex"),]
  st_cortex_exp <- as.data.frame(t(st_exp_raw[1:523, ]))
  
  sc_Seurat <- CreateSeuratObject(counts = sc_exp_raw, min.features = 10)
  sc_Seurat@meta.data$celltype <- sc_meta_raw$celltype
  sc_Seurat %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30) %>%
    Seurat::FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.3, verbose = FALSE) -> sc_Seurat
  
  st_Seurat <- CreateSeuratObject(counts = st_cortex_exp)
  st_Seurat %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30) %>%
    Seurat::FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.3, verbose = FALSE) -> st_Seurat
  
  # Annotation ST's celltype using scRNA-seq data by "FindTransferAnchors.R" Function in Seurat package
  anchors <- FindTransferAnchors(reference = sc_Seurat, query = st_Seurat)
  predictions.assay <- TransferData(anchorset = anchors, refdata = sc_Seurat$celltype, prediction.assay = TRUE, weight.reduction = st_Seurat[["pca"]], dims = 1:50)
  st_Seurat[["predictions"]] <- predictions.assay
  st_Seurat_ct <- st_Seurat[["predictions"]]@data
  for (i in 1:dim(st_Seurat_ct)[2]){
    tmp <- which(st_Seurat_ct[, i] == st_Seurat_ct['max', i], arr.ind = TRUE)[[1]]
    st_cortex_meta$celltype[i] <- rownames(st_Seurat_ct)[tmp]
  }
  
  # Select common genes in scRNA-seq and ST data
  gene_com <- intersect(rownames(sc_exp_raw), rownames(st_cortex_exp))
  if (length(gene_com) == 0) stop('The number of identical genes cannot be 0!')
  
  sc_exp <- as.data.frame(sc_exp_raw[gene_com, ])
  st_exp <- as.data.frame(st_cortex_exp[gene_com, ])
  sc_meta <- sc_meta_raw
  st_meta <- st_cortex_meta
  
  Data_Preprocessed <- list(st_exp, st_meta, sc_exp, sc_meta)
  return(Data_Preprocessed)
}




# Simulating spot-like ST data ####
#' @st_exp: DataFrame of preprocessed st data with 10000 genes in 523 cells.
#' @st_meta: DataFrame of preprocessed st metadata. 
#' @side_length: The average size of each spot.
#' @Return: A list of Simulated ST data, spot_names, spot_gc_Sum, spot_xy, spot_ct_prop

Data_Simulating <- function(st_exp, st_meta, side_length){
  
  if (!is.data.frame(st_exp)) stop('ERROR: st_exp must be a DataFrame!')
  if (!is.data.frame(st_meta)) stop('ERROR: st_meta must be a DataFrame!')
  if (!is.numeric(side_length)) stop('ERROR: st_meta must be a Numeric!')
  if (!('X' %in% colnames(st_meta))) stop('ERROR: the colnames of st_meta must have X !')
  if (!('Y' %in% colnames(st_meta))) stop('ERROR: the colnames of st_meta must have Y !')
  if (!('celltype' %in% colnames(st_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  if (!('cellname' %in% colnames(st_meta))) st_meta$cellname <- rownames(st_meta)
  
  print('Start data Simulating to construct Simulated data I.')
  
  simu_xy <- st_meta[, c('cellname', 'X', 'Y')]
  
  # The boundary and the scope of all cells
  boundary <- data.frame(x_min = min(simu_xy$X), x_max = max(simu_xy$X), 
                         y_min = min(simu_xy$Y), y_max = max(simu_xy$Y))
  
  # Compute cell sets of each spot
  Scope_x = floor(seq(boundary$x_min, boundary$x_max, by = side_length))
  Scope_y = floor(seq(boundary$y_min, boundary$y_max, by = side_length))
  cells_sxy <- list()
  cells_sxy_list <- list()
  for (sx in Scope_x) {
    cells_sx <- simu_xy$cellname[which(simu_xy$X > sx & simu_xy$X < sx+side_length)]
    # cat('The scope of sx:', sx, 'to', sx+side_length, '\n')
    for (sy in Scope_y) {
      cells_sy <- simu_xy$cellname[which(simu_xy$Y > sy & simu_xy$Y < sy+side_length)]
      # cat('The scope of sy:', sy, 'to', sy+side_length, '\n')
      cell_sxy_tmp <- intersect(cells_sx, cells_sy)
      cells_sxy_list[[which(Scope_y == sy)]] <- cell_sxy_tmp
    }
    cells_sxy[[which(Scope_x == sx)]] <- cells_sxy_list
  }
  
  # spot number
  spot_num = 0
  for (i in 1:length(cells_sxy)) {
    non_zero <- length(which(lengths(cells_sxy[[i]]) != 0))
    spot_num <- spot_num + non_zero
  }
  # cat("The number of spots:", spot_num,'\n') # side_length = 200, 85 spots
  
  
  # Name each spot
  spot_names <- NULL 
  spot_names_tmp <- NULL
  for (i in 1:length(cells_sxy)){
    for (j in 1:length(cells_sxy[[i]])){
      if (length(cells_sxy[[i]][[j]]) != 0){
        cells_nz <- cells_sxy[[i]][[j]]
        spot_names_tmp <- paste0('spot_', i , '_', j)
        spot_names <- rbind(spot_names, spot_names_tmp)
      }
    }
  }
  
  
  # Compute gene expression of each spot
  spot_gc_Sum <- NULL
  for (i in 1:length(cells_sxy)){
    for (j in 1:length(cells_sxy[[i]])){
      if (length(cells_sxy[[i]][[j]]) != 0){
        cells_nz <- cells_sxy[[i]][[j]]
        if (length(cells_sxy[[i]][[j]]) == 1){
          spot_sum <- st_exp[, cells_nz]
        }
        else{
          spot_sum <- rowSums(st_exp[, cells_nz])
        }
        spot_gc_Sum <- cbind(spot_gc_Sum, spot_sum)
      }
    }
  }
  colnames(spot_gc_Sum) <- spot_names
  
  
  # Compute spots coordinates of each spot
  spot_xy <- NULL
  for (i in 1:length(cells_sxy)){
    for (j in 1:length(cells_sxy[[i]])){
      if (length(cells_sxy[[i]][[j]]) != 0){
        cells_nz <- cells_sxy[[i]][[j]]
        if (length(cells_sxy[[i]][[j]]) == 1){
          spot_xy_tmp <- simu_xy[, c('X', 'Y')][cells_nz, ]
        }
        else{
          spot_xy_tmp <- colMeans(simu_xy[, c('X', 'Y')][cells_nz, ])
        }
        spot_xy <- rbind(spot_xy, spot_xy_tmp)
      }
    }
  }
  rownames(spot_xy) <- spot_names
  
  
  # Compute cell type proportion of each spot
  Celltypes <- unique(st_meta$celltype)
  spot_ct_prop <- matrix(nrow = spot_num, ncol = length(Celltypes), dimnames = list(spot_names[], Celltypes))
  for (i in 1:length(cells_sxy)) {
    for (j in 1:length(cells_sxy[[i]])) {
      if (length(cells_sxy[[i]][[j]]) != 0) {
        st <- paste0('spot_', i , '_', j)
        spot_ct_tmp <- st_meta[cells_sxy[[i]][[j]], 'celltype']
        
        if (length(cells_sxy[[i]][[j]]) == 1) {
          spot_ct_prop[st, spot_ct_tmp] <- 1
        } else {
          if (length(unique(spot_ct_tmp)) == 1){
            spot_ct_prop[st, unique(spot_ct_tmp)] <- 1
          } else {
            spot_ct_prop[st, names(table(spot_ct_tmp))] <- table(spot_ct_tmp) / length(spot_ct_tmp)
          }
        }
      }
      spot_ct_prop[st, is.na(spot_ct_prop[st,])] <- 0
    }
  }
  
  
  # Compute gene expression of cells from the same cell type in each spot
  genenum <- dim(st_exp)[1]
  simu_geneCT_exp <- list()
  spot_geneCT_exp <- matrix(nrow = genenum, ncol = length(Celltypes), 
                            dimnames = list(rownames(st_exp), Celltypes))
  
  for (i in 1:length(cells_sxy)) {
    for (j in 1:length(cells_sxy[[i]])) {
      if (length(cells_sxy[[i]][[j]]) != 0) {
        st <- paste0('spot_', i , '_', j)
        spot_geneCT_exp <- matrix(nrow = genenum, ncol = length(Celltypes), 
                                  dimnames = list(rownames(st_exp), Celltypes))
        
        if (length(cells_sxy[[i]][[j]]) == 1) { 
          cells_nz_ct <- st_meta[rownames(st_meta) == cells_sxy[[i]][[j]], 'celltype'] 
          spot_geneCT_exp[, cells_nz_ct] <- st_exp[, cells_sxy[[i]][[j]]]
        } else { 
          cells_nz_ct <- c()
          for (cs in 1:length(cells_sxy[[i]][[j]])) {
            cells_nz_ct[cs] <- st_meta[rownames(st_meta) == cells_sxy[[i]][[j]][cs], 'celltype'] 
          }
          spot_geneCT_exp_tmp <- st_exp[, cells_sxy[[i]][[j]]]
          if (length(unique(cells_nz_ct)) == 1) { 
            spot_geneCT_exp[, unique(cells_nz_ct)] <- rowSums(spot_geneCT_exp_tmp)
          } else {
            spot_geneCT_exp_tmp_ct <- rbind(spot_geneCT_exp_tmp, cells_nz_ct)
            rep_ct <- list()
            for (k in 1:length(colnames(spot_geneCT_exp))) {
              rep_ct[[k]] <- as.data.frame(spot_geneCT_exp_tmp_ct[ , which(spot_geneCT_exp_tmp_ct[genenum+1, ] == 
                                                                             colnames(spot_geneCT_exp)[k])])
              if (length(rep_ct[[k]]) == 1 ){
                spot_geneCT_exp[, colnames(spot_geneCT_exp)[k]] <- as.matrix(unlist(rep_ct[[k]]))[1:genenum,]
              } else {
                spot_geneCT_exp[, colnames(spot_geneCT_exp)[k]] <- rowSums(matrix(as.numeric(unlist(rep_ct[[k]][1:genenum, ])), 
                                                                                  nrow = genenum))
              }
            }
          }
        }
        spot_geneCT_exp[is.na(spot_geneCT_exp)] <- 0
        simu_geneCT_exp[[st]] <- spot_geneCT_exp
      }
    }
  }
  for (st in 1:length(simu_geneCT_exp)) {
    if (typeof(simu_geneCT_exp[[st]]) == 'character'){
      simu_geneCT_exp[[st]] <- matrix(as.numeric(simu_geneCT_exp[[st]]), ncol = length(Celltypes), 
                                      dimnames = list(rownames(st_exp), Celltypes))
    }
  }
  
  simulated_exp <- spot_gc_Sum
  simulated_xy <- spot_xy
  ground_truth <- spot_ct_prop
  #cat("The number of spots:", spot_num,'\n') 
  cat("The simulated ST data has", dim(simulated_exp)[1], 'genes in', dim(simulated_exp)[2], 'simulated spots', '\n')
  
  Data_Simulated <- list(simulated_exp, simulated_xy, ground_truth, simu_geneCT_exp)
  return(Data_Simulated)
}





# Normalization scRNA-seq and simulated st data ####
#' @data: DataFrame of gene expression data with genes by cells/spots
#' @method: Three Normalization methods ('Log_Normalization', "Zeroone_Normalization", "Sumto1_Normalization")
#' @return: Normalized gene expression matrix with genes by cells/spots

Data_Normalization <- function(Expression_data, method = c('Log_Normalization', "Zeroone_Normalization", "Sumto1_Normalization")){
  
  if (!is.data.frame(Expression_data)) stop('ERROR: Expression_data must be a DataFrame!')
  if (!(method %in% c('Log_Normalization', "Zeroone_Normalization", "Sumto1_Normalization"))) 
    stop('Please select Nomalization method from "Log_Normalization", "Zeroone_Normalization", "Sumto1_Normalization"')
  
  data <- Expression_data
  
  if (method == 'Log_Normalization'){
    print('Start Log Normalizing.')
    base = 2
    col_names <- colnames(data)
    row_names <- rownames(data)
    data <- as.matrix(data)
    if (min(data) < 0) {
      stop("Log transformation could not process negative data")
    }
    if (base == "e") {
      base = exp(1)
    } else {
      base = as.numeric(base)
    }
    data_LogNorm <- log(1+data, base = base)
    colnames(data_LogNorm) <- col_names
    rownames(data_LogNorm) <- row_names
    data_LogNorm <- as.data.frame(data_LogNorm)
    return(data_LogNorm)
  }
  
  if (method == 'Zeroone_Normalization'){
    print('Start Zero to one Normalizing.')
    col_names <- colnames(data)
    row_names <- rownames(data)
    data <- as.matrix(data)
    data_Zero1 <- matrix(nrow = length(row_names), ncol = length(col_names), dimnames = list(row_names, col_names))
    for (i in 1:length(row_names)) {
      data_Zero1[i, ] <- (data[i, ] - min(data[i, ])) / (max(data[i, ]) - min(data[i, ]))
    }
    data_Zero1 <- as.data.frame(data_Zero1)
    return(data_Zero1)
  }
  
  if (method == "Sumto1_Normalization"){
    print('Start Sum to one Normalization.')
    gene_names <- rownames(data)
    cell_names <- colnames(data)
    data <- as.matrix(data)
    data_sum_to_1 <- NULL
    for (g in 1:length(gene_names)) {
      data_sum1_tmp <- data[g,]/sum(data[g, ])
      data_sum_to_1 <- rbind(data_sum_to_1, data_sum1_tmp)
      cat("The number of gene g:", g,'\n')
    }
    colnames(data_sum_to_1) <- cell_names
    rownames(data_sum_to_1) <- gene_names
    data_sum_to_1 <- as.data.frame(data_sum_to_1)
    return(data_sum_to_1)
  }
}





Find_CTSgenes <- function(sc_exp, sc_meta){
  
  if (!is.data.frame(sc_exp)) stop('ERROR: sc_exp must be a DataFrame!')
  if (!is.data.frame(sc_meta)) stop('ERROR: sc_meta must be a DataFrame!')
  if (!('celltype' %in% colnames(sc_meta))) stop('ERROR: the colnames of st_meta must have celltype !')
  
  print('Start Find Cell-type-specific genes.')
  # FindAllMarkers() function in Seurat package
  sc_Seu <- CreateSeuratObject(counts = sc_exp)
  sc_Seu@meta.data$celltype <- sc_meta$celltype
  sc_Seu %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunUMAP(dims = 1:30) %>%
    Seurat::FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.3, verbose = FALSE) -> sc_Seu
  Seurat::Idents(object = sc_Seu) <- sc_Seu@meta.data$celltype
  sc_Seu_markers <- Seurat::FindAllMarkers(sc_Seu, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
  
  # findMarkers_one_vs_all() function in Giotto package
  instrs = createGiottoInstructions(python_path = 'D:/anaconda/envs/Giotto-env/python.exe')
  sc_Gio <- createGiottoObject(raw_exprs = sc_exp, instructions = instrs)
  sc_Gio <- normalizeGiotto(gobject = sc_Gio)
  sc_Gio@cell_metadata$leiden_clus <- as.character(sc_meta$celltype)
  sc_Gio_markers = findMarkers_one_vs_all(gobject = sc_Gio, method = 'gini', logFC = 0.5, expression_values = 'normalized', cluster_column = 'leiden_clus')
  
  # Select shared genes as cell-type-specific genes
  celltypes <- unique(sc_meta$celltype)
  sc_Seu_genes <- list()
  sc_Gio_genes <- list()
  cts_genes_tmp <- list()
  cts_genes_num <- 0
  for (ct in celltypes) {
    #print(ct)
    sc_Seu_genes[[ct]] <- sc_Seu_markers$gene[which(sc_Seu_markers$cluster == ct)]
    sc_Gio_genes[[ct]] <- sc_Gio_markers$genes[which(sc_Gio_markers$cluster == ct)]
    cts_genes_tmp[[ct]] <- intersect(sc_Seu_genes[[ct]], sc_Gio_genes[[ct]])
    cts_genes_num <- cts_genes_num + length(cts_genes_tmp[[ct]])
  }
  
  cts_genes_df <- matrix(nrow = cts_genes_num, ncol = 2, dimnames = list(seq_len(cts_genes_num), c('gene', 'celltype')))
  qishi = 1
  for (ct in celltypes) {
    #print(ct)
    gene_num <- length(cts_genes_tmp[[ct]])
    ctrep <- rep(ct, times = gene_num)
    cts_genes_df[qishi: (qishi + gene_num - 1), "celltype"] <- ctrep
    cts_genes_df[qishi: (qishi + gene_num - 1), "gene"] <- cts_genes_tmp[[ct]]
    qishi = qishi + gene_num
  }
  cts_genes <- cts_genes_df
  return(cts_genes)
}







# Construction of Simulated data I ####

data_path <- '../Data Adjustment_Datasets/Simulated data I/'

## Read raw data and preprocessing ####
SC_exp_raw <- read.csv(paste0(data_path, 'raw data/scRNA_raw_exp.csv'), header = T, row.names = 1)
SC_meta_raw <- read.csv(paste0(data_path, 'raw data/scRNA_raw_metadata.csv'), header = T, row.names = 1)
ST_exp_raw <- read.csv(paste0(data_path, 'raw data/ST_raw_exp.csv'), header = T, row.names = 1)
ST_meta_raw <- read.csv(paste0(data_path, 'raw data/ST_raw_meta.csv'), header = T, row.names = 1)

Data_Preprocessed <- Data_Preprocessing(ST_exp_raw, ST_meta_raw, SC_exp_raw, SC_meta_raw)

ST_exp <- Data_Preprocessed[[1]]
ST_meta <- Data_Preprocessed[[2]]
SC_exp <- Data_Preprocessed[[3]]
SC_meta <- Data_Preprocessed[[4]]
write.csv(SC_exp, paste0(data_path, 'simulated st/SC_exp.csv'))
write.csv(SC_meta, paste0(data_path, 'simulated st/SC_meta.csv'))

rm(ST_exp_raw, ST_meta_raw, SC_exp_raw, SC_meta_raw)


## Simulate ST data ####
ST_Simulated <- Data_Simulating(ST_exp, ST_meta, side_length = 200)

ST_Simulated_exp <- as.data.frame(ST_Simulated[[1]])
ST_Simulated_xy <- as.data.frame(ST_Simulated[[2]])
Ground_Truth <- as.data.frame(ST_Simulated[[3]])
write.csv(ST_Simulated_exp, paste0(data_path, 'simulated st/ST_Simulated_exp.csv'))
write.csv(ST_Simulated_xy, paste0(data_path, 'simulated st/ST_Simulated_xy.csv'))
write.csv(Ground_Truth, paste0(data_path, 'simulated st/Ground_Truth.csv'))

rm(ST_exp, ST_meta, Data_Preprocessed)
rm(ST_Simulated)

## Select cell-type-specific genes ####
SC_exp <- read.csv(paste0(data_path, 'simulated st/SC_exp.csv'), header = T, row.names = 1, check.names = F)
SC_meta <- read.csv(paste0(data_path, 'simulated st/SC_meta.csv'), header = T, row.names = 1, check.names = F)

CTS_genes <- Find_CTSgenes(SC_exp, SC_meta)
write.csv(CTS_genes, paste0(data_path, 'simulated st/CTS_genes.csv'))







