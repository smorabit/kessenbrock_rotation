import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from copy import deepcopy
import pickle
import sys

class velocity_helper():

	def __init__(self):
		pass

	def barcode_celltype_match(self, cell_type_file, vlm):
		 #get cell barcode : identity from cell_id_file
	    cell_id_dict = {l.split()[0].split("_")[1] : l.split()[1] for l in open(cell_type_file, 'r').readlines()}

	    #get barcodes from loom file:
	    vlm_barcodes = [b for b in vlm.ca["CellID"]]
	    
	    barcode_bool = np.array([True if barcode in cell_id_dict.keys() else False for barcode in vlm_barcodes])
	    return  np.array([cell_id_dict[barcode] for barcode in vlm_barcodes if barcode in cell_id_dict.keys()]), barcode_bool

	def processing_pipeline(self, vlm, cell_type_file=False, outputdir=False, name="velocity"):

		# remove cells with extremely low unspliced detection
		vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

		if cell_type_file != False:

			#filter based on cell_type file:
			vlm.ca["ClusterName"], cell_type_filter = self.barcode_celltype_match(
			    cell_type_file,
			    vlm
			)
			vlm.ca['CellID'] = vlm.ca['CellID'][cell_type_filter]
			try:
				vlm.filter_cells(bool_array = np.array(cell_type_filter))
			except:
				pass

		#set clusters to cell type, defined previously
		vlm.set_clusters(vlm.ca["ClusterName"])

		vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
		vlm.filter_genes(by_detection_levels=True)
		plt.savefig("{}/{}_gene_filter.pdf".format(outputdir, name))

		#normalize data by total molecule count
		vlm._normalize_S(relative_size=vlm.S.sum(0), target_size=vlm.S.sum(0).mean())
		vlm._normalize_U(relative_size=vlm.U.sum(0), target_size=vlm.U.sum(0).mean())

		# run pca and k-nearest-neighbors
		vlm.perform_PCA()
		vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)

		#fit gammas:
		vlm.fit_gammas(limit_gamma=False, fit_offset=False)

		#calculate velocity:
		vlm.predict_U()
		vlm.calculate_velocity()
		vlm.calculate_shift(assumption="constant_velocity")
		# ind4.calculate_shift(assumption="constant_unspliced", delta_t=10)
		vlm.extrapolate_cell_at_t(delta_t=1.)

		#tsne embedding
		tsne = TSNE()
		vlm.ts = tsne.fit_transform(vlm.pcs[:, :25])

		vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
		vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

		return vlm

	def plot_average_velocity_tsne(self, vlm, steps=(40,40), dpi=100):
		vlm.calculate_grid_arrows(smooth=0.8, steps=steps, n_neighbors=300, quiver_scale=0.10, minlength=0.35, min_mass=1, scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True})
		vlm.plot_grid_arrows(quiver_scale=quiver_scale,
                    scatter_kwargs_dict=scatter_kwargs_dict, min_mass=1, angles='xy', scale_units='xy',
                    headaxislength=2.75, headlength=5, headwidth=4.8, minlength=minlength,
                    plot_random=False, scale_type='absolute')

	def plot_tsne(self, vlm):
		vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:,1], c=ind4.colorandum, s=7.5)
		for i in list(set(vlm.ca["ClusterName"])):
		    ts_m = np.median(vlm.ts[vlm.ca["ClusterName"] == i, :], 0)
		    plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["ClusterName"] == i][0]),
		             fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
		plt.axis("off");

	def plot_pca(self, vlm)
		vcy.scatter_viz(vlm.pcs[:,0], vlm.pcs[:,1], c=ind4.colorandum, s=7.5)
		for i in list(set(vlm.ca["ClusterName"])):
		    ts_m = np.median(vlm.pcs[vlm.ca["ClusterName"] == i, :], 0)
		    plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["ClusterName"] == i][0]),
		             fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
		plt.axis("off");

	def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
	    new_cmap = colors.LinearSegmentedColormap.from_list(
	        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
	        cmap(np.linspace(minval, maxval, n)))
	    return new_cmap

	def scrublet_predictions(self, vlm, input_dir, doublet_rate = 0.06):
		import scrublet as scr
		import scipy.io
		print('Loading counts matrix {}/matrix.mtx'.format(input_dir), file=sys.stderr)
		counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
		print("Loading barcodes {}/barcodes.tsv".format(input_dir), file=sys.stderr)
		barcodes = np.array(scr.load_genes(input_dir + 'barcodes.tsv', delimiter='t', column=0))

		print("Initializing scrublet object", file=sys.stderr)
		scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=doublet_rate) #whole counts matrix
		print("Computing doublet predictions", file=sys.stderr)
		doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
		                                                          min_cells=3, 
		                                                          min_gene_variability_pctl=85, 
		                                                          n_prin_comps=30)
		#collapse barcodes, scores, and predictions into a dict
		doublet_dict = {barcode: [doublet_scores[i], predicted_doublets[i]] for i, barcode in enumerate(barcodes)}

		#add doublet score and doublet prediction as column attributes:
		vlm.ca["doublet_scores"] = np.array([doublet_dict[barcode][0] for barcode in vlm.ca['CellID']])
		vlm.ca["doublet_predictions"] = np.array([doublet_dict[barcode][1] for barcode in vlm.ca['CellID']])
		return vlm

	def filter_doublets(self, vlm):
		vlm.filter_cells(bool_array = np.array([not bool(doublet_prediction) for doublet_prediction in vlm.ca["doublet_predictions"]]))

	def load_marker_genes(self, marker_gene_file):
		df = pd.read_table(marker_gene_file)
		df.set_index("gene", inplace=True)
		return df

vh = velocity_helper()

ind4 = vcy.VelocytoLoom("ind4/velocyto_out/possorted_genome_bam_5AUZJ.loom")
ind4.ca['CellID'] = np.array([b.split(":")[1] for b in ind4.ca["CellID"]])

vh.processing_pipeline(vlm=ind4, cell_type_file="/home/smudge/Documents/kessenbrock_lab/RNA_velocity/ind4_celltypes.tsv")
