import sys,argparse
import os,glob
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=18
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks")
# from scipy.interpolate import interpn
# from scipy.stats import gaussian_kde



peak_files = ['UTX_peaks','UTXFEB_peaks','UTX_islands','UTXFEB_islands']
indirs = ['f2_combined_data','f3_combined_data_RPKM']

k4me1_log2fc_col='H3K4me1_WT_over_H3K4me1_Vector_log2FC'
k4me1_log2avg_col = 'H3K4me1_WT_over_H3K4me1_Vector_log2AVG'
hm_log2fc_col='H3K27ac_WT_over_H3K27ac_Vector_log2FC'
hm_log2avg_col = 'H3K27ac_WT_over_H3K27ac_Vector_log2AVG'
fc_thres = [2,1.5]
log2avg_thre = 0

for indir in indirs:
    outdir='{}_k4em1_k27ac_increased'.format(indir)
    os.makedirs(outdir,exist_ok=True)
    for peak_file in peak_files[:]:
        df = pd.read_csv('{}/combined_DiffNormReads_on_{}.csv'.format(indir,peak_file),index_col=0)
        # keep those regions with increased k4me3 in WT vs. Vector
        for fc_thre in fc_thres:
            df = df[(df[k4me1_log2fc_col]> np.log2(fc_thre)) & (df[k4me1_log2avg_col]>log2avg_thre)] 
            df_tmp = df[(df[hm_log2fc_col]> np.log2(fc_thre)) & (df[hm_log2avg_col]>log2avg_thre)] 
            df_tmp.to_csv('{}/fcthre_{}_combined_DiffNormReads_on_{}.csv'.format(outdir,fc_thre,peak_file))
    


