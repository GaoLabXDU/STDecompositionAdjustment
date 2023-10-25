
# The decomposition of DSTG methods.

# we downloaded DSTG from https://github.com/Su-informatics-lab/DSTG.

# reference:
# [27] Song Q, Su J. DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Brief Bioinform. 2021;22(5):bbaa414. doi:10.1093/bib/bbaa414.

sc_exp_path='../Datasets/Simulated data I/simulated st/SC_exp.rds'
sc_meta_path='../Datasets/Simulated data I/simulated st/SC_meta.rds'
st_exp_path='../Datasets/Simulated data I/simulated st/ST_exp.rds'

out_dir='DSTG_res/'
madir out_dir

Rscript convert_data.R sc_exp_path st_exp_path sc_meta_path
python train.py