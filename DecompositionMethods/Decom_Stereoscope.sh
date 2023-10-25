
# The decomposition of stereoscope methods.

# we downloaded stereoscope from https://github.com/almaan/stereoscope

# reference:
# [23] Andersson A, Bergenstråhle J, Asp M, Bergenstråhle L, Jurek A, Fernández Navarro J, et al. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Commun Biol. 2020;3(1):565. doi:10.1038/s42003-020-01247-y

sc_exp_path='../Datasets/Simulated data I/simulated st/SC_exp.tsv'
sc_meta_path='../Datasets/Simulated data I/simulated st/SC_meta.tsv'
st_exp_path='../Datasets/Simulated data I/simulated st/ST_exp.tsv'

out_dir='Stereoscope_res/'

stereoscope run --sc_cnt sc_exp_path --sc_labels sc_meta_path -sce 10000  \
-o out_dir -n 5000 --st_cnt st_exp_path -ste 10000 --gpu -stb 100 -scb 100
