import velocyto as vcy
import numpy as np
import pickle
from sklearn.manifold import TSNE


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

# brca = {
#     "ind1": False,
#     "ind2": True,
#     "ind3": True,
#     "ind4": True,
#     "ind9": False,
#     "ind10": False
# }
# cell_types = ["Basal", "Luminal_1", "Luminal_2"]
# brca_spliced = {cell: {acc: {"spliced": 0., "unspliced": 0.} for acc in ind.ra['Accession']} for cell in cell_types}
# normal_spliced = {cell: {acc: {"spliced": 0., "unspliced": 0.} for acc in ind.ra['Accession']} for cell in cell_types}
#
# for key in brca.keys():
#     print(key)
#     subset_bool = np.array(["{}.".format(key) in b for b in ind.ca['CellID']])
#     subset_ids = ind.ca["CellID"][subset_bool]
#     subset_clusters = ind.ca["ClusterName"][subset_bool]
#     spliced = ind.S.T[subset_bool].T
#     unspliced = ind.U.T[subset_bool].T
#     for cell in cell_types:
#         print(cell)
#         cell_bool = np.array([cell == cluster for cluster in subset_clusters])
#         cell_spliced = spliced.T[cell_bool].T
#         cell_unspliced = unspliced.T[cell_bool].T
#         for i, acc in enumerate(ind.ra['Accession']):
#             if brca[key] == True:
#                 brca_spliced[cell][acc]["unspliced"] = np.mean(cell_unspliced[i])
#                 brca_spliced[cell][acc]["spliced"] = np.mean(cell_spliced[i])
#             else:
#                 normal_spliced[cell][acc]["unspliced"] = np.mean(cell_unspliced[i])
#                 normal_spliced[cell][acc]["spliced"] = np.mean(cell_spliced[i])
#     print()

pickle.dump(ind, open("/pub/smorabit/velo/merged2.pickle", 'wb'))
# pickle.dump(brca_spliced, open("/pub/smorabit/velo/brca_spliced.pickle", 'wb'))
# pickle.dump(normal_spliced, open("/pub/smorabit/velo/normal_spliced.pickle", 'wb'))
