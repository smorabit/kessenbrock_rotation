import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from copy import deepcopy
import pickle
import sys
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from sklearn.neighbors import NearestNeighbors
import igraph
from numpy_groupies import aggregate, aggregate_np


class VelocityHelper():

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
		print("Filtering cells with low unspliced detection:", file=sys.stderr)
		vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))
		
		print("Incorporating cell types from {}".format(cell_type_file), file=sys.stderr)
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

		print("Filtering genes by detection level", file=sys.stderr)
		vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
		vlm.filter_genes(by_detection_levels=True)
		plt.savefig("{}/{}_gene_filter.pdf".format(outputdir, name))

		#normalize data by total molecule count
		print("Normalizing spliced and unspliced matrices by molecule counts", file=sys.stderr)
		vlm._normalize_S(relative_size=vlm.S.sum(0), target_size=vlm.S.sum(0).mean())
		vlm._normalize_U(relative_size=vlm.U.sum(0), target_size=vlm.U.sum(0).mean())

		# run pca and k-nearest-neighbors
		print("Running PCA", file=sys.stderr)
		vlm.perform_PCA()

		print("Running knn imputation", file=sys.stderr)
		vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)

		#fit gammas:
		print("Fitting gammas", file=sys.stderr)
		vlm.fit_gammas(limit_gamma=False, fit_offset=False)

		#calculate velocity:
		print("Calculating rna velocity", file=sys.stderr)
		vlm.predict_U()
		vlm.calculate_velocity()
		vlm.calculate_shift(assumption="constant_velocity")
		# ind4.calculate_shift(assumption="constant_unspliced", delta_t=10)
		vlm.extrapolate_cell_at_t(delta_t=1.)

		#tsne embedding
		print("Computing TSNE embedding", file=sys.stderr)
		tsne = TSNE()
		vlm.ts = tsne.fit_transform(vlm.pcs[:, :25])

		print("Estimating transition probabilities", file=sys.stderr)
		vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
							 n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
		vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

		return vlm

	def plot_raw_velocity_tsne(self, vlm, quiver_scale=25):
		plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
					c="0.8", alpha=0.2, s=10, edgecolor="")

		ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
		plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
					c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

		quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25, width=0.00045,edgecolors="k", color=vlm.colorandum[ix_choice], alpha=1)
		plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
				   vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
				   scale=quiver_scale, **quiver_kwargs)

	def plot_average_velocity_tsne(self, vlm, quiver_scale=0.1, minlength=0.35, min_mass=1, steps=(40,40), dpi=100):
		vlm.calculate_grid_arrows(smooth=0.8, steps=steps, n_neighbors=300, quiver_scale=quiver_scale, minlength=minlength, min_mass=min_mass, scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True})
		vlm.plot_grid_arrows(quiver_scale=quiver_scale,
					scatter_kwargs_dict=scatter_kwargs_dict, min_mass=min_mass, angles='xy', scale_units='xy',
					headaxislength=2.75, headlength=5, headwidth=4.8, minlength=minlength,
					plot_random=False, scale_type='absolute')

	def plot_tsne(self, vlm):
		vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:,1], c=vlm.colorandum, s=7.5)
		for i in list(set(vlm.ca["ClusterName"])):
			ts_m = np.median(vlm.ts[vlm.ca["ClusterName"] == i, :], 0)
			plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["ClusterName"] == i][0]),
					 fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
		plt.axis("off");

	def plot_doublets_tsne(self, vlm):
		vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:,1], c=vlm.ca['doublet_predictions'], s=7.5, cmap=self.truncate_colormap(plt.get_cmap('binary'),0.2))
		for i in list(set(vlm.ca["ClusterName"])):
			ts_m = np.median(vlm.ts[vlm.ca["ClusterName"] == i, :], 0)
			plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["ClusterName"] == i][0]),
					 fontsize=13, bbox={"facecolor":"w", "alpha":0.6})

	def plot_pca(self, vlm):
		vcy.scatter_viz(vlm.pcs[:,0], vlm.pcs[:,1], c=vlm.colorandum, s=7.5)
		for i in list(set(vlm.ca["ClusterName"])):
			ts_m = np.median(vlm.pcs[vlm.ca["ClusterName"] == i, :], 0)
			plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["ClusterName"] == i][0]),
					 fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
		plt.axis("off");

	def truncate_colormap(self, cmap, minval=0.0, maxval=1.0, n=100):
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

	def despline(self):
		ax1 = plt.gca()
		# Hide the right and top spines
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		# Only show ticks on the left and bottom spines
		ax1.yaxis.set_ticks_position('left')
		ax1.xaxis.set_ticks_position('bottom')

	def minimal_xticks(self, start, end):
		end_ = np.around(end, -int(np.log10(end))+1)
		xlims = np.linspace(start, end_, 5)
		xlims_tx = [""]*len(xlims)
		xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
		plt.xticks(xlims, xlims_tx)

		
	def minimal_yticks(self, start, end):
		end_ = np.around(end, -int(np.log10(end))+1)
		ylims = np.linspace(start, end_, 5)
		ylims_tx = [""]*len(ylims)
		ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
		plt.yticks(ylims, ylims_tx)

	def marker_gene_plot(self, vlm, cluster, gene_df, bad_gene_list=[], cols=2, figsize=(16.5,7.5), num_genes=6, dpi=120):

		gene_list = []
		for gene in list(gene_df.loc[gene_df['cluster'] == cluster].index):
			if gene in vlm.ra["Gene"] and gene not in bad_gene_list:
				gene_list.append(gene)
			if len(gene_list) == num_genes:
				break

		if num_genes % cols != 0:
			rows = int(num_genes/cols) + 1
		else:
			rows = int(num_genes/cols)
			
		plt.figure(None, figsize, dpi=dpi)
		gs = plt.GridSpec(rows,cols*3)
		for i, gn in enumerate(gene_list):
			ax = plt.subplot(gs[i*3])
			try:
				ix=np.where(vlm.ra["Gene"] == gn)[0][0]
			except:
				continue
			vcy.scatter_viz(vlm.Sx_sz[ix,:], vlm.Ux_sz[ix,:], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
			plt.title(gn)
			xnew = np.linspace(0,vlm.Sx[ix,:].max())
			plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
			plt.ylim(0, np.max(vlm.Ux_sz[ix,:])*1.02)
			plt.xlim(0, np.max(vlm.Sx_sz[ix,:])*1.02)
			self.minimal_yticks(0, np.max(vlm.Ux_sz[ix,:])*1.02)
			self.minimal_xticks(0, np.max(vlm.Sx_sz[ix,:])*1.02)
			self.despline()

			vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i*3+1], s=3, rasterized=True)
			vlm.plot_expression_as_color(gene_name=gn, gs=gs[i*3+2], s=3, rasterized=True)
		plt.tight_layout()

	def locate_source_sink(self, vlm, figsize=(14,7), dpi=100):
		steps = 100, 100
		grs = []
		for dim_i in range(vlm.embedding.shape[1]):
			m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
			m = m - 0.025 * np.abs(M - m)
			M = M + 0.025 * np.abs(M - m)
			gr = np.linspace(m, M, steps[dim_i])
			grs.append(gr)

		meshes_tuple = np.meshgrid(*grs)
		gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

		nn = NearestNeighbors()
		nn.fit(vlm.embedding)
		dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

		diag_step_dist = np.sqrt((meshes_tuple[0][0,0] - meshes_tuple[0][0,1])**2 + (meshes_tuple[1][0,0] - meshes_tuple[1][1,0])**2)
		min_dist = diag_step_dist / 2
		ixs = ixs[dist < min_dist]
		gridpoints_coordinates = gridpoints_coordinates[dist.flat[:]<min_dist,:]
		dist = dist[dist < min_dist]

		ixs = np.unique(ixs)

		vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='forward', cells_ixs=ixs)
		vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=2500)

		diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
		diffused_n /= np.percentile(diffused_n, 97)
		diffused_n = np.clip(diffused_n, 0, 1)

		fig, ax = plt.subplots(1,2,figsize=figsize, dpi=dpi)
		plt.subplot(121)
		vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
						c=diffused_n, alpha=0.5, s=50, lw=0.,
						edgecolor="", cmap="viridis_r", rasterized=True)
		plt.title("Sinks")
		plt.axis("off");

		vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='backwards', cells_ixs=ixs)
		vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=2500)

		diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
		diffused_n /= np.percentile(diffused_n, 97)
		diffused_n = np.clip(diffused_n, 0, 1)

		plt.subplot(122)
		vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
						c=diffused_n, alpha=0.5, s=50, lw=0.,
						edgecolor="", cmap="viridis_r", rasterized=True)
		plt.title("Sources")
		plt.axis("off");

	def array_to_rmatrix(self, X):
		nr, nc = X.shape
		xvec = robj.FloatVector(X.transpose().reshape((X.size)))
		xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
		return xr

	def principal_curve(self, X, pca=True):
		"""
		input : numpy.array
		returns:
		Result::Object
			Methods:
			projections - the matrix of the projectiond
			ixsort - the order ot the points (as in argsort)
			arclength - the lenght of the arc from the beginning to the point
		"""
		# convert array to R matrix
		xr = self.array_to_rmatrix(X)
		
		if pca:
			#perform pca
			t = robj.r.prcomp(xr)
			#determine dimensionality reduction
			usedcomp = max( sum( np.array(t[t.names.index('sdev')]) > 1.1) , 4)
			usedcomp = min([usedcomp, sum( np.array(t[t.names.index('sdev')]) > 0.25), X.shape[0]])
			Xpc = np.array(t[t.names.index('x')])[:,:usedcomp]
			# convert array to R matrix
			xr = self.array_to_rmatrix(Xpc)

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

	def pseudotime_analysis(self, vlm, figsize=(9,9), dpi=100):
		nn = NearestNeighbors(50)
		nn.fit(ind4.pcs[:,:4])
		knn_pca = nn.kneighbors_graph(mode='distance')
		knn_pca = knn_pca.tocoo()
		G = igraph.Graph(list(zip(knn_pca.row, knn_pca.col)), directed=False, edge_attrs={'weight': knn_pca.data})
		VxCl = G.community_multilevel(return_levels=False, weights="weight")
		labels = np.array(VxCl.membership)

		pc_obj = principal_curve(ind4.pcs[:,:4], False)
		pc_obj.arclength = np.max(pc_obj.arclength) - pc_obj.arclength
		labels = np.argsort(np.argsort(aggregate_np(labels, pc_obj.arclength, func=np.median)))[labels]

		manual_annotation = {str(i):[i] for i in labels}
		annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
		clusters = np.array([annotation_dict[i] for i in labels])
		colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
		vlm.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})

		k = 550
		vlm.knn_imputation(n_pca_dims=n_comps,k=k, balanced=True,
						   b_sight=np.minimum(k*8, vlm.S.shape[1]-1),
						   b_maxl=np.minimum(k*4, vlm.S.shape[1]-1))

		vlm.normalize_median()
		vlm.fit_gammas(maxmin_perc=[2,95], limit_gamma=True)

		vlm.normalize(which="imputed", size=False, log=True)
		vlm.Pcs = np.array(vlm.pcs[:,:2], order="C")

		vlm.predict_U()
		vlm.calculate_velocity()
		vlm.calculate_shift()
		vlm.extrapolate_cell_at_t(delta_t=1)

		vlm.estimate_transition_prob(hidim="Sx_sz", embed="Pcs", transform="log", psc=1,
							 n_neighbors=150, knn_random=True, sampled_fraction=1)

		vlm.calculate_grid_arrows(smooth=0.9, steps=(25, 25), n_neighbors=200)

		vlm.set_clusters(ind4.ca["ClusterName"])
		plt.figure(None,figsize=figsize, dpi=dpi)
		vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.7, "lw":0.7, "edgecolor":"0.4", "s":70, "rasterized":True},
							 min_mass=2.9, angles='xy', scale_units='xy',
							 headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.35, scale_type="absolute")
		plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="w", lw=6, zorder=1000000)
		plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="k", lw=3, zorder=2000000)
		plt.gca().invert_xaxis()
		plt.axis("off")
		plt.axis("equal");

	def marker_gene_pseudotime_plot(self, vlm, cluster, gene_df, bad_gene_list=[], cols=2, figsize=(7,6), num_genes=6, dpi=120):

		gene_list = []
		for gene in list(gene_df.loc[gene_df['cluster'] == cluster].index):
			if gene in vlm.ra["Gene"] and gene not in bad_gene_list:
				gene_list.append(gene)
			if len(gene_list) == num_genes:
				break
		
		if num_genes % cols != 0:
			rows = int(num_genes/cols) + 1
		else:
			rows = int(num_genes/cols)
			
		plt.figure(None, figsize=figsize,dpi=dpi)
		gs = plt.GridSpec(rows,cols)
		for n, gene in enumerate(gene_list):
			i = np.where(vlm.ra["Gene"] == gene)
			ax = plt.subplot(gs[n])
			plt.scatter(pc_obj.arclength[pc_obj.ixsort], vlm.Ux_sz[i, pc_obj.ixsort],
						alpha=0.7, c=np.array([0,159,193])/255, s=5, label="unspliced")
			plt.scatter(pc_obj.arclength[pc_obj.ixsort], vlm.Sx_sz[i, pc_obj.ixsort]*vlm.gammas[i],
						alpha=0.7, c=np.array([251, 172, 71])/255, s=5, label="spliced")
			m = 0 
			#np.minimum(np.min(vlm.Ux_sz[i,:]), np.min(vlm.Sx_sz[i,:]*vlm.gammas[i]))
			M = np.maximum(np.max(vlm.Ux_sz[i,:]), np.max(vlm.Sx_sz[i,:]*vlm.gammas[i]))
			plt.ylim(m - 0.07*(M-m), M + 0.07*(M-m))
			plt.ylabel(gene)
			plt.yticks([m,0.5*(m+M),M], [f"{m:.2f}", "", f"{M:.2f}"])
			p = np.min(pc_obj.arclength[pc_obj.ixsort])
			P = np.max(pc_obj.arclength[pc_obj.ixsort])
			plt.xticks(np.linspace(p,P,5), [f"{p:.0f}", "","","", f"{P:.0f}"])
			# Hide the right and top spines
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			# Only show ticks on the left and bottom spines
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('bottom')
			ax.spines['left'].set_bounds(m, M)
			ax.spines['bottom'].set_bounds(p, P)
			if n == 1:
				plt.legend()
		plt.tight_layout()

	def save_object(self, vlm, filepath):
		pickle.dump(vlm, open(filepath, 'wb'))

	def load_object(self, filepath):
		return pickle.load(open(filepath, 'rb'))
