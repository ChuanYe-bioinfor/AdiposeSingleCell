
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import scanpy as sc
import scanpy.external as sce
import os
import random

os.getcwd()
os.chdir(r'C:\path\to\Nat-Mouse')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly = True, facecolor='white')


#%% subset, preprocess
male = sc.read('male_final.h5ad')
male = male.raw.to_adata()
apc = male[male.obs.cluster == 'Adipose progenitor cell', :]

apc.raw = apc
apc.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(apc, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(apc, save = 'hvg.pdf')
apc = apc[:, apc.var.highly_variable]
sc.pp.regress_out(apc, ['pct_counts_mt'])
sc.pp.scale(apc, max_value=10)

#%% dim reduction
sc.tl.pca(apc, svd_solver='arpack', random_state=0)
sc.pl.pca_variance_ratio(apc, n_pcs=50, log=True, save = 'pca.pdf')

# batch effect
sc.pp.neighbors(apc, n_neighbors=20, n_pcs=50)
sc.tl.umap(apc)
sc.pl.umap(apc, color=['diet'], save='_apc_batch_diet.pdf')
sc.pl.umap(apc, color=['depot'], save='_apc_batch_depot.pdf')
sc.pl.umap(apc, color=['sample_id'], save='_apc_batch_sample.pdf')

# harmony
sce.pp.harmony_integrate(apc, key=['diet','depot','sample_id'], max_iter_harmony=50)
sc.pp.neighbors(apc, n_neighbors=20, n_pcs=50, use_rep='X_pca_harmony', random_state=0)
sc.tl.umap(apc, random_state=0)
sc.pl.umap(apc, color=['diet'], save='_apc_removebatch_diet.pdf')
sc.pl.umap(apc, color=['depot'], save='_apc_removebatch_depot.pdf')
sc.pl.umap(apc, color=['sample_id'], save='_apc_removebatch_sample.pdf')
sc.tl.leiden(apc, resolution=0.21)
sc.pl.umap(apc, color='leiden', save='_apc_removebatch_leiden.pdf')


#%% marker gene
sc.tl.rank_genes_groups(apc, groupby='leiden', method='wilcoxon')
markerlist_5clusters = pd.DataFrame(apc.uns['rank_genes_groups']['names']).head(50)

new_cluster_names = ['APC-1-1', 
                     'APC-2',
                     'APC-3',
                     'APC-1-2',
                     'APC-4']
apc.rename_categories('leiden', new_cluster_names)
apc.obs['cluster'] = apc.obs.leiden
apc.obs['cluster'].replace(['APC-1-1','APC-1-2'], 'APC-1', inplace=True)
sc.pl.umap(apc, color='cluster', title='', legend_loc='on data', legend_fontsize='small', legend_fontoutline=1, save='_removebatch_4cluster.pdf')

sc.tl.rank_genes_groups(apc, groupby='cluster', method='wilcoxon')
markerlist_4clusters = pd.DataFrame(apc.uns['rank_genes_groups']['names']).head(50)
for i in apc.obs.cluster.unique():
    markergene = sc.get.rank_genes_groups_df(apc, group=i)
    markergene['cluster'] = i
    markergene.to_csv('markergene_{}.csv'.format(i), index=False)

markergenes = ['Auts2', 'Bmper', 'Col4a1', 'Col4a2',
               'Creb5', 'Anxa3', 'Dpp4', 'Pi16',
               'Pid1', 'Prickle1', 'Fmo2', 'F3',
               'Top2a', 'Smc4', 'Knl1', 'Mki67']
sc.pl.dotplot(apc, var_names=markergenes, groupby='cluster', dendrogram=False, save='_apc_markergene.pdf')

new_cluster_names2 = ['Committed progenitor',
                      'Early progenitor',
                      'Regulatory progenitor',
                      'Proliferative progenitor']
apc.rename_categories('cluster', new_cluster_names2)
sc.pl.umap(apc, color='cluster', title='', legend_fontsize='small', save='_apc_removebatch_cluster_anno.pdf')

apc.write('apc_final.h5ad')
apc = apc.raw.to_adata()
apc.write('apc_forSeurat.h5ad')


# figures
sc.pl.violin(apc, ['Col4a1'], stripplot=False, groupby='cluster', save='_Col4a1.pdf')
sc.pl.violin(apc, ['Dpp4'], stripplot=False, groupby='cluster', save='_Dpp4.pdf')
sc.pl.violin(apc, ['Fmo2'], stripplot=False, groupby='cluster', save='_Fmo2.pdf')
sc.pl.violin(apc, ['Mki67'], stripplot=False, groupby='cluster', save='_Mki67.pdf')

Reds = plt.cm.get_cmap('Reds', 256)
newReds = ListedColormap(Reds(np.linspace(0.15, 1, 256)))
sc.pl.umap(apc, color='Col4a1', title='', legend_fontsize='small', color_map=newReds, save='_Col4a1.pdf')
sc.pl.umap(apc, color='Dpp4', title='', legend_fontsize='small', color_map=newReds, save='_Dpp4.pdf')
sc.pl.umap(apc, color='Fmo2', title='', legend_fontsize='small', color_map=newReds, save='_Fmo2.pdf')
sc.pl.umap(apc, color='Mki67', title='', legend_fontsize='small', color_map=newReds, save='_Mki67.pdf')


#%% Differential expressed genes of HFD vs NCD
apc = sc.read('apc_final.h5ad')
apc = apc.raw.to_adata()
apc = apc[apc.obs.leiden != 'Proliferative progenitor', :]
apc = apc[apc.obs.leiden != 'Regulatory progenitor', :]

apc.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(apc, groupby='diet', groups=['HFD'], reference='NCD', method='wilcoxon')
sc.pl.rank_genes_groups(apc, groups=['HFD'], n_genes=20)
sc.pl.rank_genes_groups_violin(apc, groups='HFD', n_genes=5)
pd.DataFrame(apc.uns['rank_genes_groups']['names']).head(25)
HFDvsNCD = sc.get.rank_genes_groups_df(apc, group='HFD')

HFDvsNCD.to_csv('HFDvsNCD_markergene.csv', index = False)









sc.pl.umap(apc, color=['Rap1a','Braf','Raf1'])







