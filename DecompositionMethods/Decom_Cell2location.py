"""
The decomposition of Cell2location methods.

We downloaded cell2location in https://pypi.org/project/cell2location/ and refered the decomposition code from [2].

References:
[17] Chen J, Liu W, Luo T, Yu Z, Jiang M, Wen J, et al. A comprehensive comparison on cell-type composition inference for spatial transcriptomics data. Brief Bioinform. 2022;23(4):bbac245. doi:10.1093/bib/bbac245.
[25] Kleshchevnikov V, Shmatko A, Dann E, Aivazidis A, King HW, Li T, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol. 2022;40(5):661–671. doi:10.1038/s41587-021-01139-4.
"""

import os
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text


out_dir = "Cell2location_res/"
output_pre = 'Simulated_data_I'
results_folder = out_dir + output_pre + "/"
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
os.makedirs(results_folder, exist_ok=True)

data_dir = "../Data Adjustment_Datasets/Simulated data I/"

scRNA_path = os.path.join(data_dir, 'simulated st/SC_exp.csv') # SC_exp data with genes by cells
st_cnt = os.path.join(data_dir, 'simulated st/ST_Simulated_exp.csv') # ST_exp data with genes by cells
scRNA_anno_path = os.path.join(data_dir, 'simulated st/SC_meta.csv')

scRna = pd.read_csv(scRNA_path, sep=',', index_col=0)
scRna = scRna.transpose() # use this step when the input is gene*cnt
anno = pd.read_csv(scRNA_anno_path, sep=",", index_col=0)

anno.index=anno.iloc[:,0].values
intersect = anno.index.intersection(scRna.index)
scRna = scRna.loc[intersect,]
anno = anno.loc[intersect,]

gene_name=pd.DataFrame(scRna.columns.values, columns=['gene'], index=scRna.columns.values)

adata_ref1 = anndata.AnnData(X=scRna, obs= anno, var=gene_name, dtype='int32')
# adata_ref1.X
# adata_ref1.var
adata_ref1.obs['Sample'] = "scRNA"

from cell2location.utils.filtering import filter_genes
selected = adata_ref1.var.index
# filter the object (not filter in our analysis)
adata_ref1 = adata_ref1[:, selected].copy()

from cell2location.models import Cell2location, RegressionModel
RegressionModel.setup_anndata(adata=adata_ref1, batch_key='Sample', labels_key=anno.columns.values[1],)
mod = RegressionModel(adata_ref1)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=5000, batch_size=512, train_size=1, lr=0.002, use_gpu=False)

# plot ELBO loss history during training, removing first 20 epochs from the plot
plt.clf()
mod.plot_history(0)
plt.savefig(results_folder+"sc_loss.png")
plt.clf()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref1 = mod.export_posterior(adata_ref1, sample_kwargs={'num_samples': 1000,
                                                             'batch_size': 512, 'use_gpu': False})


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref1.write(adata_file)
# adata_file

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref1 = sc.read_h5ad(adata_file)

run_name = f'{results_folder}/cell2location_map/'
st = pd.read_csv(st_cnt, sep=",", index_col=0)
st=st.transpose()
anno=pd.DataFrame(st.index.values,columns=['spots'],index=st.index.values)

gene_name=pd.DataFrame(st.columns.values,columns=['gene'],index=st.columns.values)
adata_vis1 = anndata.AnnData(X=st,obs= anno,var=gene_name, dtype='int32')
adata_vis1.obs['sample'] = "st"

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref1.varm.keys():
    inf_aver = adata_ref1.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                           for i in adata_ref1.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref1.var[[f'means_per_cluster_mu_fg_{i}'
                               for i in adata_ref1.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref1.uns['mod']['factor_names']
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis1.var_names, inf_aver.index)
adata_vis1 = adata_vis1[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis1, batch_key="sample")
#scvi.data.view_anndata_setup(adata_vis1)

sc.pp.filter_genes(adata_vis1, min_counts=1)
sc.pp.filter_cells(adata_vis1,min_genes=3)
# create and train the model
mod2 = cell2location.models.Cell2location(
    adata_vis1, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)
mod2.train(max_epochs= 20000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=False)

# plot ELBO loss history during training, removing first 100 epochs from the plot
plt.clf()
mod2.plot_history(1000)
plt.legend(labels=['full data training'])
plt.savefig(f'{results_folder}mod2_train.png')

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis1 = mod2.export_posterior(
    adata_vis1, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)
# Save model
mod2.save(f"{run_name}", overwrite=True)
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis1.write(adata_file)
# adata_file

plt.clf()
mod2.plot_QC()
plt.savefig(f'{results_folder}/plot_QC.png')
plt.clf()


# Obtain cell type proportions
adata=anndata.read(os.path.join(results_folder, 'cell2location_map/sp.h5ad'))
CTProp = adata.obsm['q05_cell_abundance_w_sf']
CTProp.to_csv(os.path.join(results_folder, 'cell2location_map/Cell2location_CTprop.csv'))

