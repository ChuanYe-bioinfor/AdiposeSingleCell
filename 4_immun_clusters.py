
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os

os.getcwd()
os.chdir(r'C:\path\to\Nat-Mouse')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly = True, facecolor='white')


#%% subset, preprocess
male = sc.read('male_final.h5ad')
male = male.raw.to_adata()
immun = male[(male.obs.cluster == 'Macrophage')|(male.obs.cluster == 'T cell')|(male.obs.cluster == 'B cell'),:]

immun.raw = immun
immun.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(immun, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(immun, save = 'hvg.pdf')
immun = immun[:, immun.var.highly_variable]
sc.pp.regress_out(immun, ['pct_counts_mt'])
sc.pp.scale(immun, max_value=10)


#%% dim reduction
sc.tl.pca(immun, svd_solver='arpack', random_state=0)
sc.pl.pca_variance_ratio(immun, n_pcs=50, log=True, save = 'pca.pdf')

# without harmony
sc.pp.neighbors(immun, n_neighbors=20, n_pcs=50)
sc.tl.umap(immun)
sc.pl.umap(immun, color=['diet'], save='_immun_batch_diet.pdf')
sc.pl.umap(immun, color=['depot'], save='_immun_batch_depot.pdf')
sc.pl.umap(immun, color=['sample_id'], save='_immun_batch_sample.pdf')

# harmony
sce.pp.harmony_integrate(immun, key=['diet','depot','sample_id'], max_iter_harmony=50)
sc.pp.neighbors(immun, n_neighbors=30, n_pcs=50, use_rep='X_pca_harmony', random_state=0)
sc.tl.umap(immun, random_state=0)
sc.pl.umap(immun, color=['diet'], save='_immun_removebatch_diet.pdf')
sc.pl.umap(immun, color=['depot'], save='_immun_removebatch_depot.pdf')
sc.pl.umap(immun, color=['sample_id'], save='_immun_removebatch_sample.pdf')
sc.pl.umap(immun, color='cell_type_custom', save='_immun_custom_cluster.pdf')
sc.tl.leiden(immun, resolution=0.3)
sc.pl.umap(immun, color='leiden', save='_immun_removebatch_leiden.pdf')


#%% marker gene
sc.tl.rank_genes_groups(immun, groupby='leiden', method='wilcoxon')
markerlist = pd.DataFrame(immun.uns['rank_genes_groups']['names']).head(50)

# marker selected
sc.pl.violin(immun, ['Lgals3','Cd36'], groupby='cluster')  # Macrophage-1 LAM
sc.pl.violin(immun, ['Mrc1','Rbpj'], groupby='cluster')  # Macrophage-2 PVM
sc.pl.violin(immun, ['Mki67','Kif15'], groupby='cluster')  # Macrophage-3 P-LAM
sc.pl.violin(immun, ['Ltbp1','Prg4'], groupby='cluster')  # Macrophage-4 RM
sc.pl.violin(immun, ['Arhgap26','Lyn'], groupby='cluster')  # Monocyte
sc.pl.violin(immun, ['Flt3','Clec9a'], groupby='cluster')  # DC
sc.pl.violin(immun, ['Bcl2','Prkcq'], groupby='cluster')  # 4 T
sc.pl.violin(immun, ['Bank1','Ms4a1'], groupby='cluster')  # 6 B
sc.pl.violin(immun, ['Cpa3','Kit'], groupby='cluster')

immun.obs['cluster'] = immun.obs.leiden
new_cluster_names = ['Lipid-associated macrophage',
                     'Perivascular-like macrophage',
                     'Monocyte',
                     'Proliferative LAM',
                     'T cell',
                     'Dendritic cell',
                     'B cell',
                     'Mast cell',
                     'Regulatory macrophage',]
immun.rename_categories('cluster', new_cluster_names)
new_cluster_names_ordered = ['Lipid-associated macrophage',
                             'Proliferative LAM',
                             'Perivascular-like macrophage',
                             'Regulatory macrophage',
                             'Monocyte',
                             'Dendritic cell',
                             'T cell',
                             'B cell',
                             'Mast cell']
immun.obs['cluster']=immun.obs['cluster'].cat.set_categories(new_cluster_names_ordered)
immun.uns['cluster_colors'] = ['#1f77b4', '#279e68', '#ff7f0e', '#d62728', '#aa40fc', '#8c564b','#e377c2', '#b5bd61', '#17becf']
sc.pl.umap(immun, color='cluster', title='', save='_immun_removebatch_cluster_anno.pdf')

markergenes = ['Lgals3','Cd36',
               'Mki67','Kif15',
               'Mrc1','Rbpj',
               'Ltbp1','Prg4',
               'Arhgap26','Lyn',
               'Flt3','Clec9a',
               'Bcl2','Prkcq',
               'Bank1','Ms4a1',
               'Cpa3','Kit']
sc.pl.dotplot(immun, var_names=markergenes, groupby='cluster', dendrogram=False, save='markergene.pdf')

immun.write('immun_final.h5ad')
immun = immun.raw.to_adata()
immun.write('immun_forSeurat.h5ad')



