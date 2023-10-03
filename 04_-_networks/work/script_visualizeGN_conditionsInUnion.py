"""

	script_visualizeGN_conditionsInUnion.py: 

		This script relies on script_generateSubNetwork.py. This script has found the top N sequences from each experimental
		condition and stored their connections on a set of files. Also, this script has produced the union network of all
		experimental condition's top N nodes. 

"""


# Imports: 
import numpy as np; 
import scipy
import matplotlib.pyplot as plt; 
from mpl_toolkits import mplot3d; 
import networkx as nx; 
import os, sys; 


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



# Variables to assit on loading (e.g. number of nodes, etc): 
nLoad = 1000; 
regionName = "r1"; 
dicFolder = {}; 
dicFolder["r1"] = "Reg1/"; 
dicFolder["r2"] = "Reg2/"; 
dicFolder["r3"] = "Reg3/"; 
fReadPositions = False; 

# Defining paths and names: 
dataPathBase = "./Data/"; 
subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks_nNodes-" + str(nLoad)); 
abundacesPath = os.path.join(dataPathBase, "Abundances/"); 
figPathOut = "./Pics/PlotsGN_"+str(regionName)+"_nNodes-"+str(nLoad)+'/'; 
strCall = "mkdir "+figPathOut; 
os.system(strCall); 


# Reading list of files to process: 
fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 


## Building the union network: 

# Loading edges: 
edges = np.loadtxt(os.path.join(subNetworkPath, "edges.csv"), delimiter=','); 
abundancesDict = {}; 

# Building graph: 
G = nx.Graph(); 
G.add_edges_from(edges); 
positionsFName = figPathOut + "nodePositions"; 
if (fReadPositions): 
	nodePositions = readDict(positionsFName); 
else: 
	nodePositions = nx.spring_layout(G); 
	saveDict(positionsFName, nodePositions); 

# Producing plot: 
fig = plt.figure(); 
nx.draw(G, pos=nodePositions); 
fig.savefig(figPathOut + "unionSubNetwork.png"); 


hubs_list = []

# Loop over experimental conditions that we wish to plot: 
for thisFName in allFileNames: 
	print("Processing "+thisFName); 
	colormap = []

	# Reset abundance dict: 
	for node in G.nodes: 
		abundancesDict[node] = 0

	#print(colormap)
	# Reading abundance of top nLoad nodes in experimental condition: 
	genotypeIDs_and_abundances = np.loadtxt(os.path.join(abundacesPath, thisFName+".csv"), delimiter=','); 
	thisGenotypeIDs = genotypeIDs_and_abundances[0:nLoad,0]; 
	thisAbundances = genotypeIDs_and_abundances[0:nLoad,1]; 
	for (gID, abundance) in zip(thisGenotypeIDs, thisAbundances): 
		abundancesDict[gID] = abundance; 
	
	aux = 0
	for e in abundancesDict.keys():
		if abundancesDict[e] >aux:
			aux= abundancesDict[e]
			hubs_list.append(int(e))
	hubs_list = list(set(hubs_list))		
	for node in G.nodes:	
		if int(node) in hubs_list:
			colormap.append('red')
			#print(node,'green')
		else:
			colormap.append('blue')
	# Rescaling for style -- play around with this! 
	abundancesToPlot = [5.0*np.log2(1.0+abundancesDict[node]) for node in G.nodes]; 
	# abundancesToPlot = [0.01*abundancesDict[node] for node in G.nodes]; 

	# Producing plot: 
	fig = plt.figure(); 
	nx.drawing.nx_pylab.draw_networkx_edges(G, pos=nodePositions); 
	nx.drawing.nx_pylab.draw_networkx_nodes(G,pos=nodePositions, node_color=colormap,  node_size=abundancesToPlot); 
	# nx.draw(G, pos=nodePositions, node_size=abundancesToPlot); 
	fig.savefig(figPathOut + "subNetwork_"+thisFName+".png"); 
	plt.close(fig); 


# plt.show(); 

