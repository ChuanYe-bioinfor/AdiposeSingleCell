
#%%
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

#%% import, add metadata
metadata = pd.read_csv('input/mouse_metadata.csv', index_col=[0])
adipose = sc.read_10x_mtx(
    r'D:\Nature 2022\mouse_cellranger',
    var_names='gene_symbols',
    cache=True)
adipose.var_names_make_unique()
adipose.obs = pd.merge(adipose.obs, metadata, left_index=True, right_index=True)

#%% filter: male, mt
sc.pp.filter_cells(adipose, min_genes=200)
sc.pp.filter_genes(adipose, min_cells=3)

adipose.var['mt'] = adipose.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adipose, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adipose, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.2, multi_panel=True, save = 'check_count.png')

adipose = adipose[adipose.obs.n_genes_by_counts < 6000, :]
adipose = adipose[adipose.obs.pct_counts_mt < 5, :]
male = adipose[adipose.obs.sex == 'male', :]

#%% doublet
a = []
for i in range(len(male.obs.index)):
    a.append(male.obs.diet[i] + '_' + male.obs.depot[i])
male.obs['diet_depot'] = a

sce.pp.scrublet(male, batch_key='diet_depot')
sce.pl.scrublet_score_distribution(male, figsize=(12, 4), save=('.pdf'))
male = male[male.obs.predicted_doublet == False, :]
male.write('male_rawcount.h5ad')

#%% normalize, regress_out
sc.pp.normalize_total(male, target_sum=1e4)
sc.pp.log1p(male)
male.raw = male
sc.pp.highly_variable_genes(male, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(male, save = 'hvg.pdf')
male = male[:, male.var.highly_variable]
sc.pp.regress_out(male, ['pct_counts_mt'])
sc.pp.scale(male, max_value=10)

#%% dim reduction
sc.tl.pca(male, svd_solver='arpack', random_state=0)
sc.pl.pca_variance_ratio(male, n_pcs=50, log=True, save = 'pca.pdf')

# batch effect
sc.pp.neighbors(male, n_neighbors=30, n_pcs=50)
sc.tl.umap(male)
sc.pl.umap(male, color=['diet'], save='_batch_diet.pdf')
sc.pl.umap(male, color=['depot'], save='_batch_depot.pdf')
sc.pl.umap(male, color=['sample_id'], save='_batch_sample.pdf')

# harmony
sce.pp.harmony_integrate(male, key=['diet','depot','sample_id'], max_iter_harmony=50)
sc.pp.neighbors(male, n_neighbors=30, n_pcs=50, use_rep='X_pca_harmony', random_state=0)
sc.tl.umap(male, random_state=0)
sc.pl.umap(male, color=['diet'], save='_removebatch_diet.pdf')
sc.pl.umap(male, color=['depot'], save='_removebatch_depot.pdf')
sc.pl.umap(male, color=['sample_id'], save='_removebatch_sample_id.pdf')
sc.pl.umap(male, color=['cell_type_custom'], save='_removebatch_original_celltype.pdf')
sc.tl.leiden(male, resolution=0.2)
sc.pl.umap(male, color='leiden', save='_removebatch_leiden.pdf')
male.write('male.h5ad')

#%% marker genes
sc.tl.rank_genes_groups(male, groupby='leiden', method='wilcoxon')
markerlist = pd.DataFrame(male.uns['rank_genes_groups']['names']).head(100)

new_cluster_names = ['Macrophage-1', 
                     'Adipocyte-1',
                     'Adipose progenitor cell-1',
                     'Adipose progenitor cell-2',
                     'Mesothelium',
                     'Macrophage-2',
                     'Macrophage-3',
                     'Endothelium',
                     'Epithelium-1',
                     'T cell',
                     'Epithelium-2',
                     'Lymphatic endothelium',
                     'B cell',
                     'Smooth muscle cell',
                     'Adipocyte-2',
                     'Macrophage-4',
                     ]
male.rename_categories('leiden', new_cluster_names)
male.obs['cluster'] = male.obs.leiden
for i in ['Macrophage','Adipocyte', 'Adipose progenitor cell', 'Epithelium']:
    male.obs['cluster'].replace([x for x in male.obs['cluster'].unique() if (i in x)], i, inplace=True)

new_cluster_names_ordered = ['Macrophage', 
                             'T cell',
                             'B cell',
                             'Adipose progenitor cell',
                             'Adipocyte',
                             'Epithelium',
                             'Mesothelium',
                             'Smooth muscle cell',
                             'Endothelium',
                             'Lymphatic endothelium']
male.obs['cluster']=male.obs['cluster'].cat.set_categories(new_cluster_names_ordered)
male.uns['cluster_colors'] = ['#1f77b4','#17becf','#aec7e8','#ff7f0e','#d62728',
                              '#e377c2','#b5bd61','#279e68','#aa40fc','#8c564b']
sc.pl.umap(male, color='cluster', title='', save='_removebatch_cluster_anno.pdf')

markergenes = ['Cd84', 'Adgre1',  # mac
               'Bcl2','Prkcq',  # T
               'Bank1','Ms4a1',  # B
               'Lama2', 'Pdgfra',  # apc
               'Cidec', 'Cd36',  # adipocyte
               'Ank3', 'Esr1',  # epi
               'Gpm6a', 'Upk3b',  # meso
               'Notch3', 'Postn',  # SMC
               'Pecam1', 'Flt1',  # endo
               'Reln', 'Prox1',  # LEC
               ]
sc.pl.dotplot(male, var_names=markergenes, groupby='cluster', dendrogram=False, save='markergene.pdf')

markergenes = ['Cd84',  # mac
               'Bcl2',  # T
               'Ms4a1',  # B
               'Pdgfra',  # apc
               'Cidec',  # adipocyte
               'Ank3',  # epi
               'Gpm6a',  # meso
               'Notch3',  # SMC
               'Pecam1',  # endo
               'Reln',  # LEC
               ]
sc.pl.dotplot(male, var_names=markergenes, groupby='cluster', dendrogram=False, save='markergene_onegene.pdf')

male.write('male_final.h5ad')


#%% GPCR expression in each clusters
malecluster = sc.read('male_final.h5ad')
male = sc.read('male_rawcount.h5ad')
sc.pp.normalize_total(male, target_sum=1e4)
sc.pp.log1p(male)
male.raw = male
# sc.pp.scale(male, max_value=10)
male.obs['cluster'] = malecluster.obs.cluster

male_mean = pd.DataFrame(columns=male.var_names, index=male.obs['cluster'].unique())                                                                                                 
for i in male.obs.cluster.unique(): 
    male_mean.loc[i] = male[male.obs['cluster']==i,:].X.mean(axis=0)
male_mean = male_mean.T

gpcr = pd.read_csv('input/gpcr_Nat_mouse.csv')
male_mean = male_mean[male_mean.index.isin(gpcr.symbol_mouse)]
male_mean = male_mean.sort_values(by='Adipocyte', ascending=False)
male_mean.to_csv('GPCR_expression_allClusters.csv')






