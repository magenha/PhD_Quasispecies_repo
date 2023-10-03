"""

	script_generateSubNetwork.py:

		This script generates a file containing only edges included in a given file of the whole subnetwork. By now, it
		generates the subnetwork for a single file.

"""

# Imports:
import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits import mplot3d;
import networkx as nx;
import os, sys;



# Variables to assit on loading (e.g. number of nodes, etc):
#nLoad = 100;
nLoad = -1;

# Defining paths and name variables:
dataPathBase = "./Data/";
grandNetworkPathIn = dataPathBase + "GN/";
abundacesPath = os.path.join(dataPathBase, "Abundances/");

# Path to store subnetworks:
if (nLoad<0):
	subNetworkPathOut = os.path.join(grandNetworkPathIn, "SubNetworks");
else:
	subNetworkPathOut = os.path.join(grandNetworkPathIn, "SubNetworks_nNodes-" + str(nLoad));
# Might be necessary to restart it:
strCall = "mkdir " + subNetworkPathOut;
os.system(strCall);

# Files to read the edges of the grand network:
fileNameGN = os.path.join(grandNetworkPathIn, "edges.csv");

# A list used to store all genotype IDs and retrieve only the subnetwork:
allGenotypeIDs = [];

# Looping over all filenames that we wish to process.
# From the grand network, we need to read edges that 'start and end' only on nodes within these processed files:
fIn = open(os.path.join(dataPathBase, "fileNames.csv"), 'r');
allFileNames = fIn.read().splitlines();
fIn.close();

for thisFileName in allFileNames:

	# File to store the subnetwork:
	fileNameSubNetwork = os.path.join(subNetworkPathOut, "edges_"+thisFileName+".csv");
	print("Processing: " + fileNameSubNetwork);
	strCall = "rm " + fileNameSubNetwork; # If the network existed, it is erased.
	os.system(strCall);

	# Reading genotype IDs from the abundances folders:
	genotypeIDs_and_abundances = np.loadtxt(os.path.join(abundacesPath, thisFileName+".csv"), delimiter=',')

	# genotypeIDs = genotypeIDs_and_abundances[0:nLoad,0].astype(int);
	# allGenotypeIDs += list(genotypeIDs);

	# # ACHT! Iker update. There were 2 bugs here. Now are corrected.
	if genotypeIDs_and_abundances.ndim == 1:
		genotypeIDs = []
		try:
			genotypeID = int(genotypeIDs_and_abundances[0])
			genotypeIDs += [genotypeID]
		except IndexError: # File is emtpy
			print('File '+fileNameSubNetwork+' was empty.')
			continue
	else:
		nLoad_ = len(genotypeIDs_and_abundances)
		genotypeIDs = genotypeIDs_and_abundances[0:nLoad_,0].astype(int)


	allGenotypeIDs += list(genotypeIDs)


	# First we store all edges that start with some ID from this condition in a file:
	dummyFileStart = "./startWith.csv";
	for gID in genotypeIDs:
		strToGrep = "'^" + str(gID) + ", '";
		strCall = "grep " + strToGrep + ' ' + fileNameGN + " >> " + dummyFileStart;
		os.system(strCall);

	# Then, from the edges starting with nodes from this condition, we read all edges ending in nodes from this condition:
	for gID in genotypeIDs:
		# Grep ending from starting files:
		strToGrep = "', " + str(gID) + "$'";
		strCall = "grep " + strToGrep + ' ' + dummyFileStart + ">> " + fileNameSubNetwork;
		os.system(strCall);

	# Removing the dummy file:
	strCall = "rm "+dummyFileStart;
	os.system(strCall);



# Obtaininig the Grand SubNetwork:
# 	We only need to do these if we are obtaining the top something.
# 	Otherwise, the Grand SubNetwork is the Grand Network itself.
if (nLoad>0):
	fileNameSubNetwork = os.path.join(subNetworkPathOut, "edges.csv");
	print("Obtaining grand subnetwork: " + fileNameSubNetwork);
	# If this file existed from a previous attempt, remove to avoid appending to a filed attempt:
	strCall = "rm " + fileNameSubNetwork;
	allGenotypeIDs = np.unique(allGenotypeIDs);

	# First we store all edges that start with some ID from this condition in a file:
	dummyFileStart = "./startWith.csv";
	for gID in allGenotypeIDs:
		strToGrep = "'^" + str(gID) + ", '";
		strCall = "grep " + strToGrep + ' ' + fileNameGN + " >> " + dummyFileStart;
		os.system(strCall);

	# Then, from the edges starting with nodes from this condition, we read all edges ending in nodes from this condition:
	for gID in allGenotypeIDs:
		# Grep ending from starting files:
		strToGrep = "', " + str(gID) + "$'";
		strCall = "grep " + strToGrep + ' ' + dummyFileStart + ">> " + fileNameSubNetwork;
		os.system(strCall);

	strCall = "rm "+dummyFileStart;
	os.system(strCall);
