'''
    This script runs all topologycal measures but in Networkx
'''

# Imports: 
import numpy as np; 
import scipy
import matplotlib.pyplot as plt; 
import networkx as nx; 
import os, sys; 
import pandas as pd;
import csv;
import itertools;
import pickle;

temps = [30,43]
cmap = {30:'blue', 43:'red'}

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
	temp = int(temp)
	step = int(step)
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



def calculateParametersLSR(x, y):
	"""
	Function that calculates w* and b* for a least square regression and the standard error for b*
	Inputs:alpha = 0.5, color = "green", label = "$P(k)$")
plt.plot(k_, pkBin_fit, linestyle="dashed", color = "maroon", 
		label = "fit
		x: array containing the x-values in log-scale
		y: array containing the y-values in log-scale
	Outputs:
		b: float containing the intercept term of the regression line
		w: float containing the slope of the regression line
		SEw: float containing the standard error of b*
	"""

	meanTarget = np.mean(y)
	meanFeature = np.mean(x)

	centeredTarget = y - meanTarget
	centeredFeature = x - meanFeature

	w = (centeredFeature @ centeredTarget)/(centeredFeature @ centeredFeature)

	b = meanTarget - w *  meanFeature

	# Standard error
	yHat = b + w * x
	n = x.shape[0]

	SEw = np.sqrt((1/(n - 2)) * ((np.sum(np.power((y - yHat), 2)))/(np.sum(np.power((x - meanFeature), 2)))))



	return b, w, SEw

# Variables to assit on loading (e.g. number of nodes, etc): 
nLoad =-1; 
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
	progGNPath = os.path.join(dataPathBase, "progressive_GN");
	figPathOut = "./Topology_results/";
	ProgressivefigPathOut = "./Topology_results/Progressive_GN/";
else:
	subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks_nNodes-" + str(nLoad)); 
	figPathOut = "./Topology_results_nNodes-" + str(nLoad)+'/';


strCall = "mkdir "+figPathOut; 
os.system(strCall); 


# Reading list of files to process: 
fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 
allFileNames.sort()


################################# Grand Network Analysis ###################################################

#Activate this if you wish to analyze the properties of the Grand Network
if 0==1:
	print('Analyzing final Grand Network')
	filepath = os.path.join(dataPathBase, "GN/edges.csv");
	thisFName = 'GrandNetwork'
	edges = np.loadtxt(filepath, delimiter=',', dtype=int);
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges);
	#G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
	nodes = list(G.nodes())

	print('Analyzing neighbors')
	degree = {n: nx.degree(G,n) for n in nodes}
	neighbor_degree = {n:np.array([nx.degree(G,k) for k in nx.neighbors(G,n)]).mean() for n in nodes}
	with open(f'{figPathOut}neighbors_{thisFName}.csv', 'w') as f:
		y = [degree[i] for i in nodes]
		z = [neighbor_degree[i] for i in nodes]
		f.write('id,node,degree,average_neighbor_degree\n')
		for i,n in enumerate(nodes):
			f.write(f"{i},{n},{degree[n]},{neighbor_degree[n]}\n")

	print('Analyzing Eigenvector Centrality')
	eige = nx.eigenvector_centrality(G, max_iter=1000,tol=1e-08)
		
	with open(f'{figPathOut}eigenvector_{thisFName}.csv', 'w') as f:
		y = [eige[i] for i in nodes]
		f.write('id,node,eigenvector_centrality,degree\n')
		for i,n in enumerate(nodes):
			f.write(f"{i},{n},{eige[n]},{degree[n]}\n")

	print('Analyzing Clusstering coefficients')
	lcc = nx.clustering(G)
	with open(f'{figPathOut}clustering_{thisFName}.csv', 'w') as f:
		y = [lcc[i] for i in nodes]
		f.write('id,node,LCC,degree\n')
		for i,n in enumerate(nodes):
			f.write(f"{i},{n},{lcc[n]},{degree[n]}\n")

	#Activate this if you wish to compute the betweenness centrality of the GN
	if 0==1:
		print('Analyzing Betweenness centrality')
		betw = nx.betweenness_centrality(G)
		with open(f'{figPathOut}betweenness_{thisFName}.csv', 'w') as f:
			f.write('id,node,betweenness_centrality\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{y[i]}\n")



################################# Subnetworks analysis #######################################

if 0==1:
	print('Analyzing subnetworks')
	for thisFName in allFileNames:
		t,step = analyze_fname(thisFName);
		print(f'Analyzing {thisFName}')
		#Loading edges
		filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
		edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
		#Build graph
		G = nx.Graph(); 
		G.add_edges_from(edges);
		#G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
		
		nodes = list(G.nodes())

		print('Analyzing neighbors')
		degree = {n: nx.degree(G,n) for n in nodes}
		neighbor_degree = {n:np.array([nx.degree(G,k) for k in nx.neighbors(G,n)]).mean() for n in nodes}
		with open(f'{figPathOut}Neighbors/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,degree,average_neighbor_degree\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{degree[n]},{neighbor_degree[n]}\n")
		
		print('Analyzing Clusstering coefficients')
		lcc = nx.clustering(G)
		with open(f'{figPathOut}Clustering Coefficients/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,LCC,degree\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{lcc[n]},{degree[n]}\n")
		
		print('Analyzing Eigenvector Centrality')

		eige = nx.eigenvector_centrality(G, max_iter=1000,tol=1e-08)
		
		with open(f'{figPathOut}Eigenvector_centrality/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,eigenvector_centrality\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{eige[n]}\n")
		
		
		print('Analyzing Betweenness centrality')
		betw = nx.betweenness_centrality(G)
		with open(f'{figPathOut}Betweenness_centrality/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,betweenness_centrality,degree\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{betw[n]},{degree[n]}\n")


################################# Progressive Subnetworks analysis #######################################
#Activate this to compute the progressive network measures  
if 0==1:
	print('Analyzing progressive Grand Networks')
	#Read all filenames
	fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
	allFileNames = fIn.read().splitlines(); 
	fIn.close(); 
	allFileNames.sort()
	subNetworkPath = progGNPath
	figPathOut = ProgressivefigPathOut

	
	for thisFName in allFileNames:
		t,step = analyze_fname(thisFName);
		print(f'Analyzing {thisFName}')
		#Loading edges
		filepath = os.path.join(subNetworkPath, thisFName);
		with open(filepath, 'rb') as f:
			edges = pickle.load(f)
		#Build graph
		G = nx.Graph(); 
		G.add_edges_from(edges);
		#G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
		
		nodes = list(G.nodes())

		print('Analyzing neighbors')
		degree = {n: nx.degree(G,n) for n in nodes}
		neighbor_degree = {n:np.array([nx.degree(G,k) for k in nx.neighbors(G,n)]).mean() for n in nodes}
		with open(f'{figPathOut}Neighbors/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,degree,average_neighbor_degree\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{degree[n]},{neighbor_degree[n]}\n")
		
		print('Analyzing Clusstering coefficients')
		lcc = nx.clustering(G)
		with open(f'{figPathOut}Clustering Coefficients/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,LCC,degree\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{lcc[n]},{degree[n]}\n")
		
		print('Analyzing Eigenvector Centrality')
		eige = nx.eigenvector_centrality(G, max_iter=1000,tol=1e-08)
		
		with open(f'{figPathOut}Eigenvector_centrality/edges_{thisFName}.csv', 'w') as f:
			f.write('id,node,eigenvector_centrality\n')
			for i,n in enumerate(nodes):
				f.write(f"{i},{n},{eige[n]}\n")
		
		if 0==1:
			print('Analyzing Betweenness centrality')
			betw = nx.betweenness_centrality(G)
			with open(f'{figPathOut}Betweenness_centrality/edges_{thisFName}.csv', 'w') as f:
				f.write('id,node,betweenness_centrality,degree\n')
				for i,n in enumerate(nodes):
					f.write(f"{i},{n},{betw[n]},{degree[n]}\n")



############ Network growth-acceleration ###########

#Activate this (only for growing networks)
if 1==1:
	print('Analyzing progressive Grand Networks')
	#Read all filenames
	fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
	allFileNames = fIn.read().splitlines(); 
	fIn.close(); 
	allFileNames.sort()
	subNetworkPath = progGNPath
	figPathOut = ProgressivefigPathOut

	results = {t:{} for t in temps}

	
	for thisFName in allFileNames:
		t,step = analyze_fname(thisFName);
		print(f'Analyzing {thisFName}')
		#Loading edges
		filepath = os.path.join(subNetworkPath, thisFName);
		with open(filepath, 'rb') as f:
			edges = pickle.load(f)
		#Build graph
		G = nx.Graph(); 
		G.add_edges_from(edges);
		#G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
		
		nodes = list(G.nodes())
		E = len(edges)
		V = len(nodes)

		results[t][step] = {'edges': E, 'nodes': V}
	

	#Plot
	plt.figure()
	for t in temps:
		x = list(results[t].keys())
		x.sort()
		y_E = [results[t][step]['edges'] for step in x]
		y_V = [results[t][step]['nodes'] for step in x]
		y_E_i = [y_E[i] - 2*y_V[i] for i in range(len(y_E))]  #Internal links
		d_y_E = [(y_E_i[i+1]-y_E_i[i])/(y_V[i+1]-y_V[i]) for i in range(len(y_E)-1)]
		b,w,E = calculateParametersLSR(np.log2(y_V[:-1]), np.log2(d_y_E))
		y_aster = 2**b * y_V[:-1]**w
		
		plt.plot(y_V[:-1], d_y_E, 'o-', c=cmap[t], label=f'{t}')
		plt.plot(y_V[:-1], y_aster, '--', c='black', label=f'fiting {round(w,3)}({int(E*1000)})')
	plt.legend()
	plt.ylabel('dE_i/dV')
	plt.xlabel('V')
	plt.yscale('log')
	plt.xscale('log')
	#plt.savefig('./acceleration_log.svg', dpi=300, format='svg')
	plt.savefig('./d_y_E_t_logfit.svg', dpi=300, format='svg')

	#Plot V(t) and fit
	plt.figure()
	for t in temps:
		x = list(results[t].keys())
		x.sort()
		y_E = [results[t][step]['edges'] for step in x]
		y_V = [results[t][step]['nodes'] for step in x]

		
		
		plt.plot(x, y_E, 'o-', c=cmap[t], label=f'{t}')
		
	plt.legend()
	plt.ylabel('d_y_E')
	plt.xlabel('t')
	plt.yscale('log')
	plt.xscale('log')
	#plt.savefig('./d_y_E_t_logfit.svg', dpi=300, format='svg')

	plt.show()
