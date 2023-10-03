"""

	script_visualizeGN_conditionsInUnion.py: 

		This script relies on script_visualizeGN_conditionsInUnion.py. 
        Explores each subNetwork and exports the distance matrix. 

"""


# Imports: 
import numpy as np; 
import scipy
import matplotlib.pyplot as plt; 
from mpl_toolkits import mplot3d; 
import networkx as nx; 
import os, sys; 
import pandas as pd;
import csv;
import itertools;


def saveDict(fToSave, dictToSave): 
	"""	saveDict function: 

			This function saves a dictionary into two files, one for its keys and another one for its values. This function is
			specific for the kind of dictionaries that we are working with, which contain a 2-dimensional array as values. 

			Inputs: 
				>> fToSave: Base name for the file where to save the data. 
					- From here, two files are created: $fToSave + "_keys.csv" and $fToSave + "_coord.csv". 
				>> dictToSave: A dictionary of the kind that we are dealing with in this code, which contain a 2-dimensional array
				as values. 

	"""

	fOut1 = open(fToSave + "_keys.csv", 'w'); 
	fOut2 = open(fToSave + "_coord.csv", 'w'); 
	for key in dictToSave.keys(): 
		fOut1.write(str(key)+'\n'); 
		fOut2.write(str(dictToSave[key][0])+', '+str(dictToSave[key][1])+'\n'); 
	fOut1.close(); 
	fOut2.close(); 

	return; 

def readDict(fToRead): 
	"""	readDict function: 

			This function reads a 2-D node layout that has been stored as a dictionary from the above function. 

			Inputs: 
				>> fToRead: Base name for the file from where to read the data. 
					- From here, two files are created: $fToSave + "_keys.csv" and $fToSave + "_coord.csv". 

			Returns: 
				<< newDict: A dictionary containing the name of nodes and their position for a 2D layout. 

	"""

	newDict = {}; 
	fIn1 = open(fToRead + "_keys.csv", 'r'); 
	fIn2 = open(fToRead + "_coord.csv", 'r'); 

	keys = fIn1.read().splitlines(); 
	values = np.loadtxt(fIn2, delimiter=',', converters=float); 
	for (iKey, key) in enumerate(keys): 
		newDict[int(key)] = np.array(values[iKey, :]); 

	return newDict; 

def analyze_fname(a):
	'''
	Function that looks in file_name the values of
	temperature and step
	
	Returns:  
		- temperature (Celsius degrees)
		- step
	'''
	if a[0] == 'c':
		step = int(a.split('-')[1][1:])
		t=43
	elif a[0] == '3':
		step = int(a.split('-')[2])
		t=30
	return t, step


# Variables to assit on loading (e.g. number of nodes, etc): 
nLoad = -1; 
regionName = "r1"; 
dicFolder = {}; 
dicFolder["r1"] = "Reg1/"; 
dicFolder["r2"] = "Reg2/"; 
dicFolder["r3"] = "Reg3/"; 
fReadPositions = False; 

# Defining paths and names: 
dataPathBase = "./Data/"; 

if nLoad<0:
	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks");
else:

	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks_nNodes-" + str(nLoad)); 
abundacesPath = os.path.join(dataPathBase, "Abundances/"); 
figPathOut = "./Output/Distance_matrixes_nNodes-"+str(nLoad)+'/'; 
strCall = "mkdir "+figPathOut; 
os.system(strCall); 


# Reading list of files to process: 
fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 

sCCs = {30:[], 43:[]}
sCCs_N = {30:[], 43:[]}
for thisFName in allFileNames: 
	print("Processing "+thisFName); 
    #Loading edges
	filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
	edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges); 
	t,step = analyze_fname(thisFName)

	CCs= [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
	sCCs[t].append(CCs)
	sCCs_N[t].append(CCs/np.array(CCs).sum())
sCCs = {t:list(itertools.chain(*sCCs[t])) for t in sCCs.keys()}
sCCs_N = {t:list(itertools.chain(*sCCs_N[t])) for t in sCCs_N.keys()}
plt.figure()
for t in sCCs.keys():
	bins = np.histogram_bin_edges(sCCs[t], bins='auto')
	counts,  bin2= np.histogram(sCCs[t], bins=bins)
	#plt.plot(bins[1:], counts, label=f't={t}')
	plt.scatter(bins[1:], counts, label=f't={t}')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('counts')
plt.xlabel('size of CCs')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig(f'./Output/sCCs.png')
plt.figure()
for t in sCCs.keys():
	bins = np.histogram_bin_edges(sCCs_N[t], bins='auto')
	counts,  bin2= np.histogram(sCCs_N[t], bins=bins)
	#plt.plot(bins[1:], counts, label=f't={t}')
	plt.scatter(bins[1:], counts, label=f't={t}')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('counts')
plt.xlabel('size of CCs / N')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig(f'./Output/sCCs_N.png')

