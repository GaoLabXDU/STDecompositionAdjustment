# STDecompositionAdjustment
Adjustment of scRNA-seq data to improve cell-type decomposition of spatial transcriptomics

## Description
Most sequencing-based spatial transcriptomics technologies do not achieve single-cell resolution where each captured location (spot) may contain a mixture of cells from heterogeneous cell types, and several cell-type decomposition methods have been proposed to estimate cell type proportions of each spot by integrating with scRNA-seq data. However, these existing methods ignore the distribution difference between scRNA-seq and spatial transcriptomics data, resulting in biased cell-type-specific genes derived from scRNA-seq data for decomposition. To address this issue, we develop an instance-based transfer learning framework to adjust scRNA-seq data by spatial transcriptomics data to correctly match cell-type-specific gene expression. We evaluate the effect of raw and adjusted scRNA-seq data on cell-type decomposition by eight leading decomposition methods using both simulated and real datasets. Experimental results show that data adjustment can effectively reduce distribution difference and improve the accuracy of decomposition, thus enabling for a more precise depiction on spatial organization of cell types. In summary, we highlight the significance of data adjustment in integrative analysis of scRNA-seq with spatial transcriptomics data, providing practical guidance for improved cell-type decomposition.


## Framework and Implement
The framework adjusted scRNA-seq data to correctly match the cell-type-specific gene expression and estimate the cell type composition of each spot. 
We constructed two simulated datasets and colloected four paired real datasets in our framework and conducted experiments on eight decomposition methods with raw / adjusted scRNA-seq and ST data as input. 
Our framework includes three steps:
1) we employed KMM method to adjust scRNA-seq data by ST data. Based on the cell-type-specific genes and cell types of scRNA-seq data, "Adjust_sc" function in KMM_Adjustment.py generated the adjusted scRNA-seq data.
2) we input the raw scRNA-seq data / adjusted scRNA-seq data and ST data, respectively, into eight decomposition methods (DecompositionMethods) (SPOTlight, SpatialDWLS, RCTD, CARD, STRIDE, stereoscope, cell2location, DSTG). 
* “DecompositionMethods”：
   * codes for eight decomposition methods.
3) we compared the Raw results and Adjusted results on some metrics.


## Related link
* SPOTlight (v0.1.7) https://github.com/MarcElosua/SPOTlight
* SpatialDWLS ("Giotto" v1.1.1) https://github.com/RubD/Giotto
* RCTD ("spacexr" v2.2.1) https://github.com/dmcable/spacexr
* CARD (v1.0) https://github.com/YingMa0107/CARD
* STRIDE (v0.0.2) https://github.com/wanglabtongji/STRIDE
* stereoscope https://github.com/almaan/stereoscope
* cell2location (v0.1.3) https://pypi.org/project/cell2location/
* DSTG https://github.com/Su-informatics-lab/DSTG


## Datasets description
 Download links for the original data used to generate the figures and results in the paper are listed in Datasets Summary in S1 Table. PDAC datasets and the scRNA-seq of MOB dataset were downloaded from NCBI GEO database (accession: GSM3036911, GSM4100723, GSE111672, GSE121891). Human heart dataset was downloaded from https://data.mendeley.com/datasets/mbvhhf8m62/2. The ISS data of Human heart came from https://github.com/Moldia/in_situ_seq. The ST data of MOB dataset were downloaded from https://www.spatialresearch.org/resources-published-datasets/doi-10-1126science-aaf2403/.  
