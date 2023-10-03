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

results = []
results_N = []
for thisFName in allFileNames: 
	print("Processing "+thisFName); 
    #Loading edges
	filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
	edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges); 

	CCs = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]

	results.append([analyze_fname(thisFName)[0],
					analyze_fname(thisFName)[1],
					len(CCs)]
	)

	results_N.append([analyze_fname(thisFName)[0],
					analyze_fname(thisFName)[1],
					len(CCs)/np.array(CCs).sum()]
	)

#Save nCCs

df = pd.DataFrame(results, columns=['temp', 'step', 'nCCs'])
df = df.sort_values(by=['temp', 'step'])

temps = list(set(df['temp'].tolist()))
plt.figure()
for t in temps:
	df_x = df[df['temp']==t]
	x = np.array(df_x['step'].tolist())
	y = np.array(df_x['nCCs'].tolist())
	plt.plot(x,y, label=f't={t}')
	plt.scatter(x,y)
	plt.xlabel('step')
	plt.ylabel('nCCs')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(f'./Output/nCCs.png')

#Save nCCs / N

df = pd.DataFrame(results_N, columns=['temp', 'step', 'nCCs'])
df = df.sort_values(by=['temp', 'step'])

temps = list(set(df['temp'].tolist()))
plt.figure()
for t in temps:
	df_x = df[df['temp']==t]
	x = np.array(df_x['step'].tolist())
	y = np.array(df_x['nCCs'].tolist())
	plt.plot(x,y, label=f't={t}')
	plt.scatter(x,y)
	plt.xlabel('step')
	plt.ylabel('nCCs/N')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig(f'./Output/nCCs_N.png')



