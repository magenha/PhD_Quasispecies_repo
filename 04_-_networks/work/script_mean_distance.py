"""

	script_visualizeGN_conditionsInUnion.py: 

		This script relies on script_visualizeGN_conditionsInUnion.py. 
        Explores each subNetwork and exports the distance matrix. 

"""


# Imports: 
import numpy as np; 
import math
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
 
abundacesPath = os.path.join(dataPathBase, "Abundances/"); 
if nLoad>0:
	figPathOut = "./Pics/PlotsGN_r1_nNodes-"+str(nLoad)+'/d_hist'; 
	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks_nNodes-" + str(nLoad));
else:
	figPathOut = "./Pics/PlotsGN_r1/d_hist"
	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks");
strCall = "mkdir "+figPathOut; 
os.system(strCall); 


# Reading list of files to process: 
fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 


default_cols = ['t', 'step', 'D']
df = pd.DataFrame([], columns=default_cols)

# Loop over experimental data: 
abundancesDict = {}
for thisFName in allFileNames: 
	results_hist = {}
	print("Processing "+thisFName); 
    #Loading edges
	filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
	edges = np.loadtxt(filepath, delimiter=','); 
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges); 

	#Choose the biggest component
	G = G.subgraph(max(nx.connected_components(G), key=len))

	nodes_list = list(G.nodes)
	N = len(nodes_list)
	dist_matrix = np.zeros((N,N))

	for i,n in enumerate(nodes_list):
		for j,m in enumerate(nodes_list[(i+1)%N:]):
			#find the min distance path between n and m
			d = nx.shortest_path_length(G, source=n, target=m)
			
			try:
				results_hist[d] +=1
			except:
				results_hist[d] =1
			#dist_matrix[i][j] = d
			#dist_matrix[j][i] = d
	#np.savetxt(figPathOut+'dist_matrix'+thisFName+'.txt', dist_matrix, delimiter=',')

	'''
	#For every node, extract abundance
	for node in G.nodes: 
		abundancesDict[node] = 0; 
	# Reading abundance of top nLoad nodes in experimental condition: 
	genotypeIDs_and_abundances = np.loadtxt(os.path.join(abundacesPath, thisFName+".csv"), delimiter=',', dtype=int, converters=float); 
	thisGenotypeIDs = genotypeIDs_and_abundances[0:nLoad,0]; 
	thisAbundances = genotypeIDs_and_abundances[0:nLoad,1]; 
	for (gID, abundance) in zip(thisGenotypeIDs, thisAbundances): 
		abundancesDict[gID] = abundance;
	
	n = np.array([abundancesDict[i] for i in nodes_list])
	#Compute  D
	D = (n.dot(dist_matrix.dot(np.transpose(n)))*N**2) / ((np.sum(n)**2)*N*(N-1))

	#Export result
	
	t, step = analyze_fname(thisFName)

	results.append([t, step, D])
	'''

	#Create the histogram

	#Sort the dict and normlize counts
	results_hist_sorted = {}
	norm = 0
	sum = 0
	for key in sorted(results_hist):
		results_hist_sorted[key] = results_hist[key]
		norm += results_hist_sorted[key] 
		sum += results_hist_sorted[key]*key
	for k in results_hist_sorted.keys():
		results_hist_sorted[k] = results_hist_sorted[k] /norm 
	mean_d = sum/norm
	

	t,step  = analyze_fname(thisFName)

	plt.figure()
	x = list(results_hist_sorted.keys())
	y = list(results_hist_sorted.values())
	plt.yscale('log')
	plt.plot(x,y, 'o-', label=f'temp={t} step={step}')
	#Add mean distance as a grid line
	plt.axvline(x=mean_d, ymin=10**-6, ymax=1, linestyle="--", c='red', label=f'mean={round(mean_d,1)}')
	
	plt.xlim(0.5, 15)
	plt.ylim(10**-6,1.0)
	plt.ylabel(f'p(d)')
	plt.xlabel(f'distance')
	plt.legend(loc='best')
	plt.savefig(f'{figPathOut}/{thisFName}.svg', format='svg')
'''
#plot results
df = pd.DataFrame(results, columns=['t', 'step', 'D'])
temps = list(set(df['t'].tolist()))
plt.figure()
for t in temps:
	#Extract the evolution of D
	x = np.array(df[df['t']==t]['step'].tolist())
	y = np.array(df[df['t']==t]['D'].tolist())

	plt.plot(x,y, label=f't={t}')
plt.xlabel('step')
plt.ylabel('D')
plt.legend(loc='best')
plt.savefig(figPathOut+'D.png')
plt.show()

with open(figPathOut+'D', 'w') as f:
	writer = csv.writer(f)
	for row in results:
		writer.writerow(row)
	

'''

