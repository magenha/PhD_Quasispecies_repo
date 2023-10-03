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
from sklearn.linear_model import LinearRegression


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

k = {30:{}, 43: {}}
k_bc = {30:{}, 43:{}}

selected_steps = [2,10,30,60]
temps = [30,43]
k_distr = {30:[], 43:[]}
k_distr_bc = {30:[], 43:[]}

k_distr_selected = {t:{s:[] for s in selected_steps} for t in temps}
k_distr_bc_selected = {t:{s:[] for s in selected_steps} for t in temps}

for thisFName in allFileNames: 
	print("Processing "+thisFName); 
    #Loading edges
	filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
	edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
	#Build graph
	G = nx.Graph(); 
	G.add_edges_from(edges); 
	t,step = analyze_fname(thisFName)

	#Get the mean k
	degrees = np.array([val for (node, val) in G.degree()])
	k_distr[t].append(degrees)
	k[t][step] = [degrees.mean(), degrees.std()]
	if step in selected_steps:
		k_distr_selected[t][step] = degrees

	G = G.subgraph(max(nx.connected_components(G), key=len))
	degrees = np.array([val for (node, val) in G.degree()])
	k_distr_bc[t].append(degrees)
	k_bc[t][step] = [degrees.mean(), degrees.std()]
	if step in selected_steps:
		k_distr_bc_selected[t][step] = degrees

cmap = {30: 'b', 43:'r'}

if 0 == 1:
	#Code for k mean on each step. whole and big-comp
	for t in k.keys():
		x = np.array(list(k[t].keys()))
		y = np.array([k[t][s][0] for s in k[t].keys()])
		#y_err = np.array([k[t][s][1] for s in k[t].keys()])

		plt.plot(x,y,'o-',label=f't={t}', c=cmap[t])
		#plt.errorbar(x,y,yerr=y_err)
		#plt.scatter(x,y)

		#x = np.array(list(k_bc[t].keys()))
		#y = np.array([k_bc[t][s][0] for s in k_bc[t].keys()])
		#plt.plot(x,y,'o--',mfc='none',label=f't={t} biggest component', c=cmap[t])
		#plt.errorbar(x,y,yerr=y_err)
		#plt.scatter(x,y, facecolors='none')

	plt.xlabel('step')
	plt.ylabel('<k>')
	plt.legend(loc='best')
	plt.tight_layout()
	plt.savefig('./Output/svg/k_mean.svg', dpi=300, format='svg')

if 1==1:
	#Code for k distrib, using all steps mixed. whole and big-comp
	for t in k_distr.keys():
		k_distr[t] = list(itertools.chain(*k_distr[t]))
		bins = np.histogram_bin_edges(k_distr[t], bins='auto')
		counts,  bin2= np.histogram(k_distr[t], bins=bins)
		#plt.plot(bins[1:], counts, label=f't={t}')
		plt.plot(bins[1:], counts*(bins[1]-bins[0])/np.sum(counts), 'o' , markersize=3,  c=cmap[t], label=f't={t}')
		k_distr_bc[t] = list(itertools.chain(*k_distr_bc[t]))
		#bins = np.histogram_bin_edges(k_distr_bc[t], bins='auto')
		#counts,  bin2= np.histogram(k_distr_bc[t], bins=bins)
		#plt.plot(bins[1:], counts, label=f't={t}')
		#plt.plot(bins[1:], counts, 'o', mfc='none', markersize=3, c=cmap[t], label=f't={t}')

	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('counts')
	plt.xlabel('k')
	plt.tight_layout()
	plt.legend(loc='best')
	plt.savefig(f'./Output/svg/k_distrib.svg', dpi=300, format='svg')
	plt.show()

#Code for visualize just certain degree distribs
#l = k_distr_bc_selected[30][selected_steps[1]]
#bins = np.histogram_bin_edges(k_distr_bc_selected[t], bins=10)
#counts,  bin2= np.histogram(k_distr_bc_selected[t], bins=bins)
#print()

'''
#selected_steps = [2,10,30,60]

smap = {2:'o', 10:'v', 30:'1', 60:'*'}
if 0==1:
	for t in k_distr.keys():
		for s in selected_steps:
			#k_distr_selected[t][s] = list(itertools.chain(*k_distr_selected[t][s]))
			bins = np.histogram_bin_edges(k_distr_selected[t][s], bins='auto')
			counts,  bin2= np.histogram(k_distr_selected[t][s], bins=bins)
			#plt.plot(bins[1:], counts, label=f't={t}')
			plt.plot(bins[1:], counts, 'o' , markersize=5,  marker=smap[s], c=cmap[t], label=f't={t} step {s}')
			#k_distr_bc_selected[t][s] = list(itertools.chain(*k_distr_bc_selected[t][s]))
			bins = np.histogram_bin_edges(k_distr_bc_selected[t][s], bins='auto')
			counts,  bin2= np.histogram(k_distr_bc_selected[t][s], bins=bins)
			#plt.plot(bins[1:], counts, label=f't={t}')
			plt.plot(bins[1:], counts, 'o', mfc='none', marker=smap[s], markersize=5, c=cmap[t], label=f't={t} step {s}')

	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('counts')
	plt.xlabel('k')
	plt.tight_layout()
	plt.legend(loc='best')
	plt.savefig(f'./Output/k_distrib_selected.png')
	plt.show()'''

'''
for step in range(len(k_distr[43])):

	k_list = k_distr[43][step]
	k_dict = {}
	for element in k_list:
		try:
			k_dict[element] +=1
		except:
			k_dict[element] =1

	x = np.log2(np.array(list(k_dict.keys())[1:]))
	y = np.log2(np.array([k_dict[e] for e in list(k_dict.keys())[1:] ]))
	#print(y)
	#x = np.reshape(x,(x.shape,1))
	#y = np.reshape(y,(y.shape,1))
	x = x.reshape(-1, 1)
	y = y.reshape(-1, 1)
	reg = LinearRegression().fit(x, y)
	print(f'step {step} score={reg.score(x, y)}, slope={reg.coef_[0]}')

	#print(reg.intercept_)
#plt.scatter(x,y)

#plt.show()'''