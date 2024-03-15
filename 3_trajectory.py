
#%%
import numpy as np
import pandas as pd
import statsmodels
import scanpy as sc
import scanpy.external as sce
import cellrank as cr
import scvelo as scv
import palantir
import os

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm

os.getcwd()
os.chdir(r'C:\path\to\Nat-Mouse')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly = True, facecolor='white', color_map = 'viridis', scanpy=True)


#%% rerun umap
apc = sc.read('apc_final.h5ad')
apc = apc.raw.to_adata()
apc = apc[apc.obs.diet == 'NCD', :]
apc = apc[apc.obs.depot == 'eWAT', :]
apc = apc[apc.obs.cluster != 'Proliferative progenitor', :]
apc = apc[apc.obs.cluster != 'Regulatory progenitor', :]

apc.raw = apc
apc.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(apc, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.tl.pca(apc, svd_solver='arpack', random_state=0)
sce.pp.harmony_integrate(apc, key=['sample_id'], max_iter_harmony=50)
sc.pp.neighbors(apc, n_neighbors=20, n_pcs=50, use_rep='X_pca_harmony', random_state=0)
sc.tl.umap(apc, random_state=0)
sc.pl.umap(apc, color=['sample_id'], save='_removebatch_sample.pdf')
sc.pl.umap(apc, color=['cluster'], title = '', legend_loc='on data', legend_fontsize='small', legend_fontoutline=1.5, save='NCD_eWAT_label.pdf')

#%% dpt
umap = pd.DataFrame(apc.obsm['X_umap'], index=apc.obs_names, columns=['x','y'])
start_cell = 'Mm_EPI_19-1_CTCAGAAGTAGAATGT'
palantir.plot.highlight_cells_on_tsne(umap, start_cell)
cell_ids = pd.Series(np.arange(len(apc.obs.index)), index=apc.obs.index)
cell_num = cell_ids['Mm_EPI_19-1_CTCAGAAGTAGAATGT']

apc.uns['iroot'] = np.flatnonzero(apc.obs['leiden'])[cell_num]
sc.tl.dpt(apc)
sc.pl.umap(apc, color=['dpt_pseudotime'], legend_fontoutline=1.5, legend_fontsize='small', legend_loc='on data', save = 'dpt.pdf')

apc.write('temp.h5ad')
apc = sc.read('temp.h5ad')


#%% sctour
# get the raw UMI counts
male = sc.read('male_rawcount.h5ad')
cells = apc.obs.index
rawapc = male[cells,:]
rawapc.var['highly_variable'] = apc.var['highly_variable']
rawapc.obs['cluster'] = apc.obs['cluster']

rawapc.write('rawapc.h5ad')
rawapc = sc.read('rawapc.h5ad')

# run sctour
rawapc = rawapc[:, rawapc.var.highly_variable]
import sctour as sct
tnode = sct.train.Trainer(rawapc, loss_mode='nb')
tnode.train()
rawapc.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
rawapc.obsm['X_TNODE'] = mix_zs
rawapc.obsm['X_VF'] = tnode.get_vector_field(rawapc.obs['ptime'].values, rawapc.obsm['X_TNODE'])

rawapc = rawapc[np.argsort(rawapc.obs['ptime'].values), :]
sc.pp.neighbors(rawapc, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(rawapc, min_dist=0.2)
sc.pl.umap(rawapc, color='cluster', size=20, legend_loc='on data', legend_fontoutline=1.5, legend_fontsize='small', save = '_latent_space_cluster.pdf')
rawapc.obs['ptime'] = sct.train.reverse_time(rawapc.obs['ptime'].values)
sc.pl.umap(rawapc, color='ptime', size=20, save = '_latent_space_ptime.pdf')

fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(18, 9))
sct.vf.plot_vector_field(rawapc, ax=axs[0], reverse=True, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', 
                         color='cluster', legend_loc='none', title='Trajectory', frameon=False, size=300, alpha=0.25, stream_linewidth=1.3)
sct.vf.plot_vector_field(rawapc, ax=axs[1], reverse=True, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', 
                         color='cluster', legend_loc='none', title='Trajectory', frameon=True, size=300, alpha=0.25, stream_linewidth=1.3)
fig.savefig('trajectory.pdf')

rawapc.write('rawapc_ptime.h5ad')


#%% cellrank_connectivity
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
scv.settings.verbosity = 3
scv.settings.set_figure_params("scanpy")
cr.settings.verbosity = 2

ck = ConnectivityKernel(apc).compute_transition_matrix()
g = GPCCA(ck)
g.compute_schur(n_components=20)

g.compute_macrostates(n_states=2, cluster_key="cluster")
plt.rcParams.update({'axes.grid': False})
g.plot_macrostates(size = 50, save='macrostate.pdf')
g.plot_coarse_T(save = 'transition_matrix.pdf')

g.set_terminal_states_from_macrostates(names=["Committed progenitor"])
g.compute_absorption_probabilities()
g.plot_absorption_probabilities(same_plot=False, size=10, basis="umap")

drivers = g.compute_lineage_drivers(lineages="Committed progenitor", return_drivers=True)
drivers = drivers.dropna(how='all')
drivers.rename(columns={'Committed progenitor corr':'corr',
                        'Committed progenitor pval':'pval',
                        'Committed progenitor qval':'FDR',
                        'Committed progenitor ci low':'CI-low',
                        'Committed progenitor ci high':'CI-high'}, inplace=True)
drivers = drivers.sort_values(by="pval", ascending=True)
pval_array = drivers['pval'].values
fdr = statsmodels.stats.multitest.multipletests(pval_array, method='fdr_bh', is_sorted=True)
drivers['FDR'] = fdr[1]
drivers = drivers.sort_values(by="corr", ascending=False)
drivers.to_csv('drivers.csv')

apc_ptime = sc.read('rawapc_ptime.h5ad')
apc.obs['ptime'] = apc_ptime.obs['ptime']

g.compute_terminal_states()
g.compute_absorption_probabilities()
model = cr.ul.models.GAM(apc)
up = list(drivers.index)[0:500]
down = list(reversed(drivers.index))[0:500]
cr.pl.heatmap(
    apc,
    model,
    genes=(up + down),
    show_absorption_probabilities=False,
    lineages="Committed progenitor",
    n_jobs=1,
    time_key="ptime",
    backend="loky",
    save = 'top500.pdf'
)

apc.obs['Committed_memberships'] = np.array(g.macrostates_memberships[:,1]).flatten()
sc.pl.umap(apc, color=['Committed_memberships'], color_map=plt.cm.Blues, save = 'macrostates_memberships.pdf')




