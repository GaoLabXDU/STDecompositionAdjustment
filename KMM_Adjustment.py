
"""
The KMM method dealt with the issue of correcting sample selection bias using unlabeled data in [1].
By estimating the probability density distribution, it reuses instances based on the weight generation rules and aims to make the distribution of the weighted source domain and target domain as close as possible in transfer learning.

Reference:
[1] Bernhard S, John P, Thomas H. Correcting Sample Selection Bias by Unlabeled Data. Advances in Neural Information Processing Systems 19: Proceedings of the 2006 Conference. MIT Press, 2007, pp.601-608.
"""

import pandas as pd
import numpy as np
import sklearn.metrics
from cvxopt import matrix, solvers
# import matplotlib.pyplot as plt
# import seaborn as sns
# import scipy.stats as stats
# from sklearn import svm
# import os
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.metrics import accuracy_score
# import argparse


def kernel(ker, X1, X2, gamma):
    K = None
    if ker == 'linear':
        if X2 is not None:
            K = sklearn.metrics.pairwise.linear_kernel(np.asarray(X1), np.asarray(X2))
        else:
            K = sklearn.metrics.pairwise.linear_kernel(np.asarray(X1))
    elif ker == 'rbf':
        if X2 is not None:
            K = sklearn.metrics.pairwise.rbf_kernel(np.asarray(X1), np.asarray(X2), gamma)
        else:
            K = sklearn.metrics.pairwise.rbf_kernel(np.asarray(X1), None, gamma)
    return K



def KMM(Xs, Xt, kernel_type, gamma=1.0, B=1.0, eps=None):
    """
    Fit source and target domain using KMM (compute the coefficients)
    :param Xs: source domain with ns samples * dim.
    :param Xt: target domain with nt samples * dim.
    :param kernel_type: kernel type, 'linear' or 'rbf'.
    :param gamma: kernel bandwidth for rbf kernel.
    :param B: bound for beta.
    :param eps: bound for sigma_beta.
    :return: probability density ratio of samples.
    """

    ns = Xs.shape[0] # sample numbers
    nt = Xt.shape[0]

    if eps == None:
        eps = B / np.sqrt(ns)

    K = kernel(kernel_type, Xs, None, gamma)
    kappa = np.sum(kernel(kernel_type, Xs, Xt, gamma) * float(ns) / float(nt), axis=1)

    K = matrix(K.astype(np.double))
    kappa = matrix(kappa.astype(np.double))
    G = matrix(np.r_[np.ones((1, ns)), -np.ones((1, ns)), np.eye(ns), -np.eye(ns)])
    h = matrix(np.r_[ns * (1 + eps), ns * (eps - 1), B * np.ones((ns,)), np.zeros((ns,))])

    sol = solvers.qp(K, -kappa, G, h)
    beta = np.array(sol['x'])
    return beta



def Adjust_sc(sc_exp, sc_meta, st_exp, cts_genes, kernel_type='rbf', gamma=None, B=2, eps=None):
    """
    Altered the weights of the corresponding cells in each cell type on the corresponding cell-type-specific genes.
    :param sc_exp: DataFrame of scRNA-seq gene expression profile with rows being genes and columns being cells.
    :param sc_meta: DataFrame of cell meta information with 'celltype'.
    :param st_exp: DataFrame of ST gene expression profile with rows being genes and columns being spots.
    :param cts_genes: DataFrame of cell-type-specific genes for each cell type.
    :param kernel_type: kernel type, 'linear' or 'rbf'.
    :param gamma: kernel bandwidth for rbf kernel.
    :param B: bound for beta.
    :param eps: bound for sigma_beta
    :return: DataFrame of adjusted scRNA-seq data with rows being genes and columns being cells of each cell type.
    """
    # Check parameters
    if not isinstance(sc_exp, pd.DataFrame):
        print("Please enter sc_exp data of pandas DataFrame type with rows being genes and columns being cells.")
    if not isinstance(sc_meta, pd.DataFrame):
        print("Please enter sc_meta data of pandas DataFrame type with cell type.")
    if not isinstance(st_exp, pd.DataFrame):
        print("Please enter st_exp data of pandas DataFrame type with rows being genes and columns being spots.")
    if not isinstance(cts_genes, pd.DataFrame):
        print("Please enter cts_genes data of pandas DataFrame type with cell type.")
    if (kernel_type != 'linear') and (kernel_type != 'rbf'):
        print("Please enter 'linear' or 'rbf' for KMM method as kernel type.")
    if not isinstance(B, float):
        print('Please enter bandwidth B of float type.')

    if 'celltype' not in sc_meta.columns:
        print("Please check the column names of sc_meta data, using 'celltype' to represent the cell types information.")
    if 'gene' not in cts_genes.columns or 'celltype' not in cts_genes.columns:
        print("Please check the column names of sc_meta data, using 'gene' and 'celltype' to represent the cell-type-specific genes and cell types.")

    print("Data Adjustment.")

    sc_expT = sc_exp.T
    st_expT = st_exp.T
    CTS = cts_genes['celltype'].unique()

    Beta_B = KMM(Xs=sc_expT, Xt=st_expT, kernel_type=kernel_type, gamma=gamma, B=B, eps=eps)
    Beta_B_pd = pd.DataFrame(Beta_B, columns=['Beta_B'])
    Beta_B_pd.index = sc_expT.index

    sc_adjust_li = []

    for ct in CTS:
        # print(ct)
        ct_cells = sc_meta.index[sc_meta['celltype'] == ct]
        ct_genes = cts_genes['gene'][cts_genes['celltype'] == ct]
        ct_Beta = Beta_B_pd.iloc[Beta_B_pd.index.isin(ct_cells.tolist())]
        ct_Beta_n = ct_Beta.copy()
        ct_sc_exp = sc_expT.loc[ct_cells, ct_genes]  # 该ct的 cells的ct_genes的基因表达

        for ce in ct_cells:
            # print(ce)
            if ct_Beta.loc[ce, 'Beta_B'] < 1:
                ct_Beta_n.loc[ce, 'Beta_B'] = ct_Beta.loc[ce, 'Beta_B'] + 1
            else:
                ct_Beta_n.loc[ce, 'Beta_B'] = ct_Beta.loc[ce, 'Beta_B']

        ct_Beta_n = ct_Beta_n.values
        ct_sc_exp_B = ct_Beta_n * ct_sc_exp
        sc_adjust_li.append(ct_sc_exp_B)

    sc_adjust_all = pd.DataFrame()
    for ls in range(0, len(sc_adjust_li)):
        diff_gen = (sc_expT.columns).difference(sc_adjust_li[ls].columns)
        same_cell = sc_adjust_li[ls].index
        sc_all_tmp = pd.concat([sc_adjust_li[ls], sc_expT.loc[same_cell, diff_gen]], axis=1)
        sc_adjust_all = pd.concat([sc_adjust_all, sc_all_tmp], axis=0)

    sc_adjust_all_Tmp = sc_adjust_all.copy()
    sc_Adjusted = sc_adjust_all_Tmp.loc[sc_expT.index.tolist(), sc_expT.columns.tolist()].T
    print('The Adjusted scRNA-seq data has : %d genes * %d cells' % (sc_Adjusted.shape))
    return sc_Adjusted







# Data Adjustment for Simulated data I
import pandas as pd
import os
import KMM_Adjustment

# Simulated data I
data_path = '../Datasets/Simulated data I/simulated st/'
SC_exp = pd.read_csv(os.path.join(data_path, 'SC_exp.csv'), sep=',', index_col=0)
SC_meta = pd.read_csv(os.path.join(data_path, 'SC_meta.csv'), sep=',', index_col=0)
ST_Simulated_exp = pd.read_csv(os.path.join(data_path, 'ST_Simulated_exp.csv'), sep=',', index_col=0)
CTS_genes = pd.read_csv(os.path.join(data_path, 'CTS_genes.csv'), sep=',', index_col=0)

Adjusted_sc = KMM_Adjustment.Adjust_sc(sc_exp=SC_exp, st_exp=ST_Simulated_exp, sc_meta=SC_meta, cts_genes = CTS_genes, kernel_type='rbf', gamma=None, B=2.0, eps=None)
Adjusted_sc.to_csv(os.path.join(data_path, 'SC_exp_Adjusted.csv'), sep=',', index=True, header=True)


