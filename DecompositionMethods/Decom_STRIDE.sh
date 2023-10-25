
# The decomposition of STRIDE methods.

# We downloaded STRIDE from https://github.com/wanglabtongji/STRIDE and follwed the instructions in https://stridespatial.readthedocs.io/en/latest/index.html#welcome-to-stride-s-documentation.

# reference:
# [26] Sun D, Liu Z, Li T, Wu Q, Wang C. STRIDE: accurately decomposing and integrating spatial transcriptomics using single-cell RNA sequencing. Nucleic Acids Res. 2022;50(7):e42. doi:10.1093/nar/gkac150

sc_exp_path='../Datasets/Simulated data I/simulated st/SC_exp.txt'
sc_meta_path='../Datasets/Simulated data I/simulated st/SC_meta.txt'
st_exp_path='../Datasets/Simulated data I/simulated st/ST_exp.txt'

out_dir='STRIDE_res/'

out_prefix=SimulatedI


STRIDE deconvolve --sc-count sc_exp_path \
--sc-celltype sc_meta_path \
--st-count st_exp_path \
--outdir out_dir --outprefix out_prefix --normalize
