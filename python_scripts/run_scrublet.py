#!/usr/bin/env python3
#
# Author: Samuel Morabito
# Kai Kessenbrock lab
# University of California, Irvine
# Updated 11/28/18
#
# refer to this link for installation instructions for Scrublet:
# https://github.com/AllonKleinLab/scrublet

import matplotlib.pyplot as plt
import numpy as np
import scrublet as scr
import scipy.io
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='raw 10X file directory for input', type=str)
parser.add_argument('-o', '--output', help='output directory', type=str, default="./")
parser.add_argument('-n', '--name', help='name of output files', type=str, default="name")
parser.add_argument('-r', '--doublet', help='expected doublet rate, default=0.06', type=float, default=0.06)
parser.add_argument('-e', '--embed', help='plot UMAP and TSNE. True or False.', type=bool, default=False)
args = parser.parse_args()

#load counts matrix, genes, barcodes
print("Loading counts matrix %s" % args.input + '/matrix.mtx', file=sys.stderr)
counts_matrix = scipy.io.mmread(args.input + '/matrix.mtx').T.tocsc()
print("Loading barcodes %s" % args.input + '/barcodes.tsv', file=sys.stderr)
barcodes = np.array(scr.load_genes(args.input + 'barcodes.tsv', delimiter='t', column=0))

#initialize scrublet object
print("Initializing scrublet object", file=sys.stderr)
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.doublet) #whole counts matrix

print("Computing doublet predictions", file=sys.stderr)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

#write scrublet output to file:
print("Writing doublet predictions to %s" % args.output + "/" + args.name + "_predicted_doublets.tsv", file=sys.stderr)
with open(args.output + "/" + args.name + "_predicted_doublets.tsv", 'w') as outfile:
	outfile.write("\t".join(["barcode", "doublet_score", "doublet_prediction"])+"\n")
	for barcode, score, prediction in zip(barcodes, doublet_scores, predicted_doublets):
		if prediction == False:
			doublet = "0"
		else:
			doublet = "1"
		outfile.write("\t".join([barcode, str(score), doublet])+"\n")

print("Plotting doublet score histogram to %s" % args.output + "/" + args.name + "_score_histogram.pdf", file=sys.stderr)
f = scrub.plot_histogram()
plt.savefig(args.output + "/" + args.name + "_score_histogram.pdf")

if args.embed == True:
	print("Running UMAP", file=sys.stderr)
	scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
	print("Plotting UMAP to %s" % args.output + "/" + args.name + "_UMAP.pdf", file=sys.stderr)
	f = scrub.plot_embedding('UMAP', order_points=True);
	plt.savefig(args.output + "/" + args.name + "_UMAP.pdf")

	print("Running TSNE", file=sys.stderr)
	scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))
	print("Plotting TSNE to %s" % args.output + "/" + args.name + "_TSNE.pdf", file=sys.stderr)
	f = scrub.plot_embedding('tSNE', order_points=True);
	plt.savefig(args.output + "/" + args.name + "_TSNE.pdf")



