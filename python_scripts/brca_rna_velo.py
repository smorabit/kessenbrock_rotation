%matplotlib inline
import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from copy import deepcopy
import pickle
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from sklearn.neighbors import NearestNeighbors
import igraph #couldn't load
from numpy_groupies import aggregate, aggregate_np #couldn't load
import loompy


################################################################################
# Step 0: Combine several loom objects                                         #
################################################################################
files = ["/pub/smorabit/velo/possorted_genome_bam_82GD0.loom", \
         "/pub/smorabit/velo/possorted_genome_bam_QKKJQ.loom", \
         "/pub/smorabit/velo/possorted_genome_bam_ZRA06.loom", \
         "/pub/smorabit/velo/possorted_genome_bam_5AUZJ.loom", \
         "/pub/smorabit/velo/possorted_genome_bam_GZ7KV.loom", \
         "/pub/smorabit/velo/possorted_genome_bam_DZN3A.loom"]
loompy.combine(files, "/pub/smorabit/velo/merged.loom", key="Accession")

################################################################################
# Step 1: Process raw data                                                     #
################################################################################

# load merged loom file:
ind = vcy.VelocytoLoom("/pub/smorabit/velo/merged.loom")
ind_name = "merged"
print_dir = "/pub/smorabit/velo/figures/"

# %%read in metadata file from Seurat
metadata = pd.read_table("/pub/smorabit/velo/Norm.BRCA.Combined.Seurat.Meta.Data.Object.txt")
metadata.set_index("barcode", inplace=True)

# rename barcodes to match seurat metadata:
prefix_conversion = {
    "possorted_genome_bam_82GD0": "ind1.epithelial.NORMAL_",
    "possorted_genome_bam_QKKJQ": "ind2.epithelial.BRCA_",
    "possorted_genome_bam_ZRA06": "ind3.epithelial.BRCA_",
    "possorted_genome_bam_5AUZJ": "ind4.epithelial.BRCA_",
    "possorted_genome_bam_GZ7KV": "ind9.epithelial.NORMAL_",
    "possorted_genome_bam_DZN3A": "ind10.epithelial.NORMAL_"
}
initial_barcodes = ind.ca['CellID']
updated_barcodes = ["{}{}".format(prefix_conversion[b.split(":")[0]], b.split(":")[1]) for b in ind.ca['CellID']]
ind.ca["CellID"] = np.array(updated_barcodes)

# how many barcodes match between merged loom object and seurat metadata in each individual?
metadata_match = metadata.loc[updated_barcodes].dropna()

#ind4 is a special snowflake so do something different for that
ind4 = open("/data/users/smorabit/ind4_celltypes.tsv")
ind4_celltypes = {"{}{}".format("ind4.epithelial.BRCA_" ,l.split()[0].split("_")[1]) : l.split()[1] for l in ind4.readlines()}
for key, val in ind4_celltypes.items():
    if val == "Basal_Myoepithelial":
        ind4_celltypes[key] = "Basal"
    elif val == "Luminal_1_1" or val == "Luminal_1_2":
        ind4_celltypes[key] = "Luminal_1"

ind4_barcodes = [b for b in updated_barcodes if b in ind4_celltypes.keys()]

# there is something wrong with individual 4 but I think I can fix it, for now I will move on since there is a kitty in my lap
# later go look on my desktop and find the notebook where I did the ind4 analysis and grab that metadata file just for ind4

# remove cells that are not found in the seurat metadata:
metadata_bool = np.array([True if (b in metadata_match.index or b in ind4_barcodes) else False for b in ind.ca["CellID"]])
ind.filter_cells(bool_array = metadata_bool)

# set cell type clusters based on the seurat metadata:
clusternames = np.array([metadata.loc[b]["Cell.Type"] if "ind4" not in b else ind4_celltypes[b] for b in ind.ca["CellID"]])
ind.ca["ClusterName"] = clusternames
ind.set_clusters(ind.ca["ClusterName"])

# %% remove cells with extremely low unspliced detection
ind.filter_cells(bool_array=ind.initial_Ucell_size > np.percentile(ind.initial_Ucell_size, 0.5))

# %%feature selection using Support Vector Regression
ind.score_detection_levels(min_expr_counts=40, min_cells_express=30)
ind.filter_genes(by_detection_levels=True)
ind.score_cv_vs_mean(3000, plot=False, max_expr_avg=35)
ind.filter_genes(by_cv_vs_mean=True) #this line doesn't work in a basic python shell

# %%normalize data by total molecule count
ind._normalize_S(relative_size=ind.S.sum(0), target_size=ind.S.sum(0).mean())
ind._normalize_U(relative_size=ind.U.sum(0), target_size=ind.U.sum(0).mean())

# %% run pca and k-nearest-neighbors
ind.perform_PCA()

# %%plot pca variance explained
# plt.plot(np.cumsum(ind.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(ind.pca.explained_variance_ratio_))>0.002))[0][0]
# plt.axvline(n_comps, c="k")
# plt.xlabel("PC")
# plt.ylabel("cumulative variance explained")
# plt.savefig(print_dir + "pca_variance_explained.pdf")

ind.knn_imputation(n_pca_dims=n_comps, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)

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
pickle.dump(ind_pseudotime, open("/pub/smorabit/velo/{}/{}_pseudotime.hdf5".format(ind_name,ind_name), 'wb'))


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

################################################################################
#                      Spliced vs unspliced                                    #
################################################################################
brca = {
    "ind1": False,
    "ind2": True,
    "ind3": True,
    # "ind4": True,
    "ind9": False,
    "ind10": False
}
brca_spliced = {acc: np.array([0.,0.]) for acc in ind.ra['Accession']}
normal_spliced = {acc: np.array([0.,0.]) for acc in ind.ra['Accession']}
for key in brca.keys():
    print(key)
    subset_bool = np.array(["{}.".format(key) in b for b in ind.ca['CellID']])
    spliced = ind.S.T[subset_bool].T
    unspliced = ind.U.T[subset_bool].T
    for i, acc in enumerate(ind.ra['Accession']):
        if brca[key] == True:
            brca_spliced[acc][0] = sum(unspliced[i])
            brca_spliced[acc][1] = sum(spliced[i])
        else:
            normal_spliced[acc] = sum(unspliced[i])
            normal_spliced[acc] = sum(spliced[i])

################################################################################
brca = {
    "ind1": False,
    "ind2": True,
    "ind3": True,
    # "ind4": True,
    "ind9": False,
    "ind10": False
}
cell_types = ["Basal", "Luminal_1", "Luminal_2"]
brca_spliced = {cell: {acc: {"spliced": 0., "unspliced": 0.} for acc in ind.ra['Accession']} for cell in cell_types}
normal_spliced = {cell: {acc: {"spliced": 0., "unspliced": 0.} for acc in ind.ra['Accession']} for cell in cell_types}

for key in brca.keys():
    print(key)
    subset_bool = np.array(["{}.".format(key) in b for b in ind.ca['CellID']])
    subset_ids = ind.ca["CellID"][subset_bool]
    subset_clusters = ind.ca["ClusterName"][subset_bool]
    spliced = ind.S.T[subset_bool].T
    unspliced = ind.U.T[subset_bool].T
    for cell in cell_types:
        print(cell)
        cell_bool = np.array([cell == cluster for cluster in subset_clusters])
        cell_spliced = spliced.T[cell_bool].T
        cell_unspliced = unspliced.T[cell_bool].T
        for i, acc in enumerate(ind.ra['Accession']):
            if brca[key] == True:
                brca_spliced[cell][acc]["unspliced"] = sum(cell_unspliced[i])
                brca_spliced[cell][acc]["spliced"] = sum(cell_spliced[i])
            else:
                normal_spliced[cell][acc]["unspliced"] = sum(cell_unspliced[i])
                normal_spliced[cell][acc]["spliced"] = sum(cell_spliced[i])
    print()

pickle.dump(ind, open("/pub/smorabit/velo/merged2.pickle", 'wb'))
pickle.dump(brca_spliced, open("/pub/smorabit/velo/brca_spliced.pickle", 'wb'))
pickle.dump(normal_spliced, open("/pub/smorabit/velo/normal_spliced.pickle", 'wb'))
