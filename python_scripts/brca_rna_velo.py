%matplotlib inline
import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from copy import deepcopy
import pickle
from velocity_helper import VelocityHelper
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
from IPython.display import display, Markdown, Latex
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from sklearn.neighbors import NearestNeighbors
import igraph
from numpy_groupies import aggregate, aggregate_np

################################################################################
#                                Load data                                     #
################################################################################

# %% Initialization
# %%load a pickled object (one that has already been processed)
ind = pickle.load(open("/home/smudge/Documents/kessenbrock_lab/RNA_velocity/ind4/ind4.hdf5", 'rb'))

# %%load raw data file (needs to be processed)
ind = vcy.VelocytoLoom("/home/smudge/Documents/kessenbrock_lab/RNA_velocity/rna_velocity_data/ind10/possorted_genome_bam_DZN3A.loom")

ind_name = "ind10"
print_dir = "/home/smudge/Documents/kessenbrock_lab/RNA_velocity/figures/{}/".format(ind_name)

################################################################################
#                            Process raw data                                  #
################################################################################

# %%read in metadata file from Seurat
metadata = pd.read_table("/home/smudge/Documents/kessenbrock_lab/RNA_velocity/Norm.BRCA.Combined.Seurat.Meta.Data.Object.txt")
metadata.set_index("barcode", inplace=True)
metadata.head()

# %%adjust barcodes to only contain nucleic acid information:
ind.ca['CellID'] = np.array([b.split(":")[1] for b in ind.ca["CellID"]])

# %% remove cells with extremely low unspliced detection
ind.filter_cells(bool_array=ind.initial_Ucell_size > np.percentile(ind.initial_Ucell_size, 0.5))

# %%filter cells based on the Seurat metadata dataframe
# %%also set cluster names based on Seurat metadata dataframe
cell_id_dict = {barcode.split("_")[1]: metadata[metadata["individual"] == ind_name]["individual.analysis.idents"][barcode] for barcode in list(metadata[metadata["individual"] == ind_name].index)}
barcode_bool = np.array([True if barcode in cell_id_dict.keys() else False for barcode in ind.ca['CellID']])
clusternames = np.array([cell_id_dict[barcode] for barcode in ind.ca['CellID'] if barcode in cell_id_dict.keys()])
ind.filter_cells(bool_array = np.array(barcode_bool))
ind.ca["ClusterName"] = clusternames

# %%filter out stromal cells if present:
stromal_filter = np.array([False if "stromal" in cell_type else True for cell_type in ind.ca["ClusterName"]])
ind.filter_cells(bool_array=stromal_filter)

# %%set clusters to cell type, defined previously. filter out cells from the unknown cluster
ind.set_clusters(ind.ca["ClusterName"])
ind.filter_cells(bool_array = np.array([True if cell_id_dict[barcode] != "Unknown" else False for barcode in ind.ca["CellID"]]))

# %%feature selection using Support Vector Regression
ind.score_detection_levels(min_expr_counts=40, min_cells_express=30)
ind.filter_genes(by_detection_levels=True)
ind.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
ind.filter_genes(by_cv_vs_mean=True)
plt.savefig(print_dir + "feature_selection.pdf")

# %%normalize data by total molecule count
ind._normalize_S(relative_size=ind.S.sum(0), target_size=ind.S.sum(0).mean())
ind._normalize_U(relative_size=ind.U.sum(0), target_size=ind.U.sum(0).mean())

# %% run pca and k-nearest-neighbors
ind.perform_PCA()
# %%ind.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=1000, b_maxl=1500, n_jobs=16)
ind.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)

# %%plot pca variance explained
plt.plot(np.cumsum(ind.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(ind.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.xlabel("PC")
plt.ylabel("cumulative variance explained")
plt.savefig(print_dir + "pca_variance_explained.pdf")


# %%calculate velocity
ind.fit_gammas() # %%there are a ton of parameters here that I could change
ind.predict_U() # %% U_predicted = gamma * S, gamma is degredation rate
ind.calculate_velocity() # %%calculates velocity as U_measured - U_predicted
ind.calculate_shift(assumption="constant_velocity") # %%compute change in gene expression for each cell
ind.extrapolate_cell_at_t(delta_t=1.) # %% extrapolates gene expression profile of each cell after a specified delta_t

# %%compute TSNE using top 50 pcs
tsne = TSNE()
ind.ts = tsne.fit_transform(ind.pcs[:, :50])

ind.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
ind.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

################################################################################
#                            Plot embeddings                                   #
################################################################################

# %%raw velocity plot
f, ax = plt.subplots(1,1, figsize=(6,6), dpi=150)
plt.subplot(111)
quiver_scale = 25
plt.scatter(ind.embedding[:, 0], ind.embedding[:, 1],
            c="0.8", alpha=0.2, s=10, edgecolor="")

ix_choice = np.random.choice(ind.embedding.shape[0], size=int(ind.embedding.shape[0]/1.), replace=False)
plt.scatter(ind.embedding[ix_choice, 0], ind.embedding[ix_choice, 1],
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25, width=0.00045,edgecolors="k", color=ind.colorandum[ix_choice], alpha=1)
plt.quiver(ind.embedding[ix_choice, 0], ind.embedding[ix_choice, 1],
           ind.delta_embedding[ix_choice, 0], ind.delta_embedding[ix_choice, 1],
           scale=quiver_scale, **quiver_kwargs)
plt.axis("off");
plt.savefig(print_dir + "raw_velocity.pdf")

# %% velocity tsne
ind.calculate_grid_arrows(smooth=0.8, steps=(30, 30), n_neighbors=300)
f, ax = plt.subplots(1,1, figsize=(10,10))
ind.plot_grid_arrows(quiver_scale=0.05,
                    scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=1, angles='xy', scale_units='xy',
                    headaxislength=2.75, headlength=5, headwidth=4.8, minlength=0.35,
                    plot_random=False, scale_type='absolute')

plt.axis("off");
plt.savefig(print_dir + "average_velocity.pdf")

# %%tsne with cluster names
f, ax = plt.subplots(1,1, figsize=(10,10))
vcy.scatter_viz(ind.ts[:,0], ind.ts[:,1], c=ind.colorandum, s=7.5)
for i in list(set(ind.ca["ClusterName"])):
    ts_m = np.median(ind.ts[ind.ca["ClusterName"] == i, :], 0)
    plt.text(ts_m[0], ts_m[1], str(ind.cluster_labels[ind.ca["ClusterName"] == i][0]),
             fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
plt.axis("off");

################################################################################
#                           Pseudotime projection                              #
################################################################################

#make a copy of ind object for pseudotime Analysis
ind_pseudotime = deepcopy(ind)

def array_to_rmatrix(X):
    nr, nc = X.shape
    xvec = robj.FloatVector(X.transpose().reshape((X.size)))
    xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
    return xr

def principal_curve(X, pca=True):
    # convert array to R matrix
    xr = array_to_rmatrix(X)

    if pca:
        #perform pca
        t = robj.r.prcomp(xr)
        #determine dimensionality reduction
        usedcomp = max( sum( np.array(t[t.names.index('sdev')]) > 1.1) , 4)
        usedcomp = min([usedcomp, sum( np.array(t[t.names.index('sdev')]) > 0.25), X.shape[0]])
        Xpc = np.array(t[t.names.index('x')])[:,:usedcomp]
        # convert array to R matrix
        xr = array_to_rmatrix(Xpc)

    #import the correct namespace
    d = {'print.me': 'print_dot_me', 'print_me': 'print_uscore_me'}
    princurve = importr("princurve", on_conflict='warn')

    #call the function
    fit1 = princurve.principal_curve(xr)

    #extract the outputs
    class Results:
        pass
    results = Results()
    results.projections = np.array( fit1[0] )
    results.ixsort = np.array( fit1[1] ) - 1 # R is 1 indexed
    results.arclength = np.array( fit1[2] )
    results.dist = np.array( fit1[3] )

    if pca:
        results.PCs = np.array(xr) #only the used components

    return results

nn = NearestNeighbors(50)
nn.fit(ind_pseudotime.pcs[:,:4])
knn_pca = nn.kneighbors_graph(mode='distance')
knn_pca = knn_pca.tocoo()
G = igraph.Graph(list(zip(knn_pca.row, knn_pca.col)), directed=False, edge_attrs={'weight': knn_pca.data})
VxCl = G.community_multilevel(return_levels=False, weights="weight")
labels = np.array(VxCl.membership)

pc_obj = principal_curve(ind_pseudotime.pcs[:,:4], False)
pc_obj.arclength = np.max(pc_obj.arclength) - pc_obj.arclength
labels = np.argsort(np.argsort(aggregate_np(labels, pc_obj.arclength, func=np.median)))[labels]

manual_annotation = {str(i):[i] for i in labels}
annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
clusters = np.array([annotation_dict[i] for i in labels])
colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
#ind_pseudotime.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})
ind_pseudotime.set_clusters(ind_pseudotime.ca["ClusterName"])

k = 550
ind_pseudotime.knn_imputation(n_pca_dims=n_comps,k=k, balanced=True,
                   b_sight=np.minimum(k*8, ind_pseudotime.S.shape[1]-1),
                   b_maxl=np.minimum(k*4, ind_pseudotime.S.shape[1]-1))

ind_pseudotime.normalize_median()
ind_pseudotime.fit_gammas(maxmin_perc=[2,95], limit_gamma=True)
ind_pseudotime.normalize(which="imputed", size=False, log=True)
ind_pseudotime.Pcs = np.array(ind_pseudotime.pcs[:,:2], order="C")
ind_pseudotime.predict_U()
ind_pseudotime.calculate_velocity()
ind_pseudotime.calculate_shift()
ind_pseudotime.extrapolate_cell_at_t(delta_t=1)
ind_pseudotime.estimate_transition_prob(hidim="Sx_sz", embed="Pcs", transform="log", psc=1,
                             n_neighbors=150, knn_random=True, sampled_fraction=1)
ind_pseudotime.calculate_grid_arrows(smooth=0.9, steps=(25, 25), n_neighbors=200)
plt.figure(None,(9,9))
ind_pseudotime.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.7, "lw":0.7, "edgecolor":"0.4", "s":70, "rasterized":True},
                     min_mass=2.9, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.35, scale_type="absolute")
plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="w", lw=6, zorder=1000000)
plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="k", lw=3, zorder=2000000)
plt.gca().invert_xaxis()
plt.axis("off")
plt.axis("equal");
plt.savefig(print_dir + "pseudotime.pdf")

#pickle processed ind_pseudotime object
pickle.dump(ind_pseudotime, open("/data1/rna_velocity_data/{}/{}_pseudotime.hdf5".format(ind_name,ind_name), 'wb'))


################################################################################
#                      Average velocity violin plots                           #
################################################################################

# %% average velocity for each cell type violin plot
mycolors = [list(thing) for thing in set(tuple(x) for x in ind.colorandum)]
mycolors = [mycolors[1], mycolors[2], mycolors[0]]

ind_splice_dict = {cl: defaultdict(list) for cl in list(set(ind.ca["ClusterName"]))}
for cluster in list(set(ind.ca["ClusterName"])):
    for i in range(ind.S.T[ind.ca['ClusterName'] == cluster].T.shape[1]):
        ind_splice_dict[cluster]["spliced"].append(np.sum(ind.S.T[ind.ca['ClusterName'] == cluster].T[:,i]))
        ind_splice_dict[cluster]["unspliced"].append(np.sum(ind.U.T[ind.ca['ClusterName'] == cluster].T[:,i]))

ratios = []
for i, cluster in enumerate(ind_splice_dict.keys()):
    ratios.append(np.divide(ind_splice_dict[cluster]['unspliced'], ind_splice_dict[cluster]['spliced']))

#compute the mean of velocities for all cells:
velocity_mean = np.array([np.mean(ind.velocity[:,i]) for i in range(ind.velocity.shape[1])])

#split this up into each cluster:
data = [velocity_mean[ind.ca["ClusterName"] == cluster] for cluster in list(set(ind.ca["ClusterName"]))]

fig, axes = plt.subplots(1,2, figsize=(6,3),dpi=150)

ax = plt.subplot(121)
violin = plt.violinplot(data, vert=False, showmeans=True);
for i, v in enumerate(violin['bodies']):
    v.set_facecolor(mycolors[i])
    v.set_edgecolor('k')
    v.set_alpha(0.7)

# change the line color from blue to black
for partname in ('cbars','cmins','cmaxes', 'cmeans'):
    vp = violin[partname]
    vp.set_edgecolor('k')
    vp.set_linewidth(1)

plt.yticks([], [])
plt.title("velocity")

#edit spines:
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
labels = [cluster for cluster in list(set(ind.ca["ClusterName"]))]
plt.yticks([i for i in range(1, len(labels)+1)], labels);
plt.tight_layout()

ax = plt.subplot(122)
violin = plt.violinplot(ratios, vert=False, showmeans=True);
for i, v in enumerate(violin['bodies']):
    v.set_facecolor(mycolors[i])
    v.set_edgecolor('k')
    v.set_alpha(0.7)

# change the line color from blue to black
for partname in ('cbars','cmins','cmaxes', 'cmeans'):
    vp = violin[partname]
    vp.set_edgecolor('k')
    vp.set_linewidth(1)

plt.yticks([], [])
plt.title("unspliced / spliced")

#edit spines:
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(print_dir + "velocity_magnitudes.pdf")
