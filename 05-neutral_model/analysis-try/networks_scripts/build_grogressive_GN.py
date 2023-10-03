'''
    This script uses the edges archives from the Data/GN folder and create progressive GrandNetworks. 
    
    The Base is to build GNs with Henry's Grand Networks too:
        -concatenating both temperatures
        -steps = 4, 10, 30, 60, ancestor(step=2), GN (all together)
'''


import numpy as np; 
import scipy
import matplotlib.pyplot as plt; 
import networkx as nx; 
import os, sys; 
import pandas as pd;
import csv;
import itertools;
import pickle


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

def redact_fname(temp, step):
	'''
        Function that gives the filename, having the temperature and step desired
    '''
	if step<10:
		step_string = '0' + str(step)
	elif step>=10:
		step_string = str(step)
	
	if temp == 30:
		x='30-1-'+step_string+'-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed_Nrem_rem'
		return x
	elif temp == 43:
		x='c43-p'+step_string+'-3-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed_Nrem_rem'
		return x
	else:
		print(f'No temperature {temp}')
		return None

# Variables to assit on loading (e.g. number of nodes, etc): 
nLoad =-1; 
regionName = "r1"; 
dicFolder = {}; 
dicFolder["r1"] = "Reg1/"; 
dicFolder["r2"] = "Reg2/"; 
dicFolder["r3"] = "Reg3/"; 
fReadPositions = False; 

temps = [30,43]
cmap = {30:'blue', 43:'red'}

# Defining paths and names: 
dataPathBase = "./Data/"; 

if nLoad<0:
	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks");
	figPathOut = "./Data/progressive_GN/";
else:

	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks_nNodes-" + str(nLoad)); 
	figPathOut = "./Data/progressive_GN"+str(nLoad)+'/'; 
abundacesPath = os.path.join(dataPathBase, "Abundances/"); 

strCall = "mkdir "+figPathOut; 
os.system(strCall); 


# Reading list of files to process: 
fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 
allFileNames.sort()



allFileNames.sort()

temp = 30

thisFName = redact_fname(temp, 2)


#Load Grand Network
print('Loading Grand Network')
filepath = os.path.join(dataPathBase, "GN/edges.csv");
thisFName = 'GrandNetwork'
edges = np.loadtxt(filepath, delimiter=',', dtype=int);
#Build graph
GN = nx.Graph(); 
GN.add_edges_from(edges);
#G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
	
nodes_GN = list(GN.nodes())


#For a given step, concatenate both temperatures



#Develop this

ab_path = dataPathBase + 'Abundances/'
temp=43
cummulative_nodes = []
index_step = 0
for thisFName in allFileNames:
	t,step = analyze_fname(thisFName)
	this_step_list = []

	filename_out = figPathOut + thisFName.replace('edges_', '')
	
	if t == temp:
		print(f'Analyzing {thisFName}')
		filename = ab_path+thisFName + '.csv'
		filename.replace('edges_', '')
		
		df = pd.read_csv(filename, header=None, names=['node', 'abundance', 'rel', 'norm'])
		nodes_step = df['node'].tolist()
		#Searching for nodes so far
		print('Searching nodes so far')
		if index_step:
			this_step_list = cummulative_nodes[index_step-1]
		elif not index_step:
			this_step_list = []
		for n in nodes_step:
			if index_step:
				if n not in cummulative_nodes[index_step-1] and n in nodes_GN:
					this_step_list.append(n)
				else:
					pass 
			elif not index_step:
				if n in nodes_GN:
					this_step_list.append(n)
				else:
					pass 
		#Find subgraph
		print('Computing subgraph')
		subgraph = GN.subgraph(this_step_list)
		print(len(nodes_step), len(list(subgraph.nodes())))
		index_step +=1
		cummulative_nodes.append(this_step_list)
		print('Saving')
		with open(figPathOut+thisFName, 'wb') as f:
			pickle.dump(list(subgraph.edges()), f)
	

		

				

          
        
#Mixes the progressive GN for 30 and 43 Celcius degrees.
#Result should be the same as Luis's