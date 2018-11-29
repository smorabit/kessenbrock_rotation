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

		print(vlm.S.shape)

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

vh = velocity_helper()

ind4 = vcy.VelocytoLoom("ind4/velocyto_out/possorted_genome_bam_5AUZJ.loom")
ind4.ca['CellID'] = np.array([b.split(":")[1] for b in ind4.ca["CellID"]])

vh.processing_pipeline(vlm=ind4, cell_type_file="/home/smudge/Documents/kessenbrock_lab/RNA_velocity/ind4_celltypes.tsv")