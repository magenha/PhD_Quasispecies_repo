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
nLoad = 100; 
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

temps = [30,43]

mean_assort = {t:{} for t in temps}
mean_assort_bc = {t:{} for t in temps}


for thisFName in allFileNames: 
	print("Processing "+thisFName); 
    #Loading edges
	filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
	edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges); 
	t,step = analyze_fname(thisFName)

	assort = nx.degree_pearson_correlation_coefficient(G)
	mean_assort[t][step] = assort

	G = G.subgraph(max(nx.connected_components(G), key=len))

	assort = nx.degree_pearson_correlation_coefficient(G)
	mean_assort_bc[t][step] = assort




cmap = {30: 'b', 43:'r'}

plt.figure()
for t in temps:
	x = list(mean_assort[t].keys())
	y = [mean_assort[t][s] for s in x]
	plt.plot(x,y,'o-',mfc='none',label=f't={t}', c=cmap[t])

plt.xlabel('step')
plt.ylabel('assortativity')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('./Output/assort')


plt.figure()
for t in temps:
	x = list(mean_assort_bc[t].keys())
	y = [mean_assort_bc[t][s] for s in x]
	plt.plot(x,y,'o-',mfc='none',label=f't={t} biggest component', c=cmap[t])

plt.xlabel('step')
plt.ylabel('assortativity')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('./Output/assort_bc')
plt.show()

