'''
	This script is intended to export a plot of the mean degree in network, 
	using different cutoff values.

'''


import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle

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
	else:
		return None
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


def calculateParametersLSR(x, y, base=2, log_transform=True):
	"""
	Function that calculates w* and b* for a least square regression and the standard error for b*
	Inputs:
		x: array containing the x-values in log-scale
		y: array containing the y-values in log-scale
	Outputs:
		b: float containing the intercept term of the regression line
		w: float containing the slope of the regression line
		SEw: float containing the standard error of b*
	"""
	if log_transform:
		#Remove the points with value -inf  
		y_aster = []
		x_aster = []
		for i,value in enumerate(y):
			if value:
				#print(value)
				y_aster.append(value)
				x_aster.append(x[i])
			else:
				print('Avoiding', i)

		
		x = np.array(x_aster)
		y = np.array(y_aster)
		if base ==2:
			x = np.log2(x)
			y = np.log2(y)
		elif base==10:
			x = np.log10(x)
			y = np.log10(y)
		else:
			x = np.log(x)
			y = np.log(y)
	x = np.array(x)
	y = np.array(y)
	'''
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
	'''
	
	m, c = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y, rcond=None)[0]

	return m,c


def alphaNewmanClauset(degrees, kmin):
	"""
	Function that calculates the exponent coefficient (alpha) for the degree distribution of a network
	Inputs: 
		degrees: array of ints containing the degrees of a network
		kmin: int indicating the minimum degree for which the power law holds
	Outputs:
		alpha: float containing the exponent of the power law
			alpha = 1 + N * (sum_{i}^{N}log((k_i)/(k_min - 0.5)))^-1
		sigma: float containing the statistical error of the estimation of alpha 
			sigma = (alpha - 1)/(srqt(N))
	"""

	# Get the number of nodes with a degree greater than or equal to kmin
	#N = degrees[degrees >= kmin].shape[0]
	degrees = np.array(degrees)
	k_min = degrees[degrees >= kmin]
	N = k_min.shape[0]

	# Perform the sum
	summ = np.sum(np.log((k_min)/(kmin - 0.5)))

	# Get alpha
	alpha = 1 + N * (1/summ)
	# Get sigma
	sigma = (alpha - 1)/np.sqrt(N)

	return alpha, sigma

#top_res_path = './Topology_results/'
top_res_path = './Topology_results/Progressive_GN/'
temps = [30, 43]
cmap = {30:'blue', 43:'red'}
symb_dict = {2:"o" , 4: "^", 10: "s", 30: "p", 60: "H"}
greek_letterz=[chr(code) for code in range(945,970)]

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

#Results of
#GN : number of nodes, number of edges, and number of nodes &edges in Giant component
#Repeat the same for each step and temp
#Save a file like this:
#GN, -, -, fraction,-, -, fraction
#t_step, -, -, fraction,-, -, fraction

#Activate this to know if we have accelerated networks.
if 0==1:
    print('#############Starting acceleration analysis#####################')
    ############ Degree #############
    n_nodes = {t:{} for t in temps}
    n_edges = {t:{} for t in temps}
    
    fig, axs = plt.subplots(1,3, figsize=(10, 5))
    #Read the graph
    allFileNames = os.listdir(progGNPath)
    allFileNames = [a for a in allFileNames if '_' in a]
    allFileNames.sort()

    for thisFileName in allFileNames:
        print(thisFileName)
        t, step = thisFileName.split('_')
        step = int(step)
        t = int(t)
        print(f'Analizing Progressive GN for t={t} and step={step}')
        thisFName = redact_fname(temp=t, step=step)
        
        with open(f'{progGNPath}/{thisFileName}', 'rb') as f:
            edges_list = pickle.load(f)
        
        G = nx.Graph()
        G.add_edges_from(edges_list)
        n_nodes[t][step] = len(list(G.nodes()))
        n_edges[t][step] = len(list(G.edges()))
	
    for t in temps:
        L = list(n_nodes[t].keys())
        L.sort()
        axs[0].plot(np.array(L), np.array([n_nodes[t][step] for step in L]), c=cmap[t], label=f'{t}ºC')
        axs[1].plot(np.array(L), np.array([n_edges[t][step] for step in L]), c=cmap[t], label=f'{t}ºC')
        axs[2].plot(np.array(L), np.array([n_edges[t][step]/n_nodes[t][step] for step in L]), c=cmap[t], label=f'{t}ºC')
    
    axs[0].set_xlabel('step')
    axs[0].set_ylabel('number of nodes')
    axs[1].set_xlabel('step')
    axs[1].set_ylabel('number of edges')
    axs[2].set_yscale('log')
    axs[2].set_xscale('log')
    axs[2].set_xlabel('step')
    axs[2].set_ylabel('number of edges/number of nodes')
    axs[0].legend(loc='best')
    axs[1].legend(loc='best')
    axs[2].legend(loc='best')
    fig.tight_layout()

    plt.savefig(f'{top_res_path}/Figures/acceleration_log.svg', dpi=300, format='svg')
    plt.show()
		    
        #G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])



#Activate this to export the evolution of mean degree and std.
if 1==1:
    print('#############Starting plot of mean degrees#####################')
    ############ Degree #############
    
    for k_min in range(1,6):
        fig, axs = plt.subplots(1,2, figsize=(10, 5))
        #Read the graph
        allFileNames = os.listdir(f'{top_res_path}Neighbors/')
        allFileNames.sort()
        desired_steps = [4, 10, 30, 60]

        #mean k vs step for both temps in a plot
        #std-dev of k vs step for both temps in a plot
        results = {t:{} for t in temps}
        
        for thisFileName in allFileNames:
            try:
                t, step = analyze_fname(thisFileName.replace('edges_', ''))
            except:
                t, step = analyze_fname(thisFileName.replace(' edges_', ''))

            df = pd.read_csv(f'{top_res_path}Neighbors/{thisFileName}')
            df = df[df['degree'] >= k_min]
            results[t][step] = {'mean': df['degree'].mean(), 'std': df['degree'].std()}
            
        
        for t in temps:
            L = list(results[t].keys())
            L.sort()

            x = np.array(L)
            y = []
            y_err = []
            for step in L:
                y.append(results[t][step]['mean'])
                y_err.append(results[t][step]['std'])
            y = np.array(y)
            y_err = np.array(y_err)
            axs[0].plot(x,y, 'o-', c=cmap[t],  label=f'{t}ºC, k_min={k_min}')
            axs[1].plot(x,y_err, 'o-', c=cmap[t],  label=f'{t}ºC')
        axs[0].set_xlabel('step')
        axs[0].set_ylabel('<k>')
        #axs[0].set_xscale('log')
        #axs[0].set_yscale('log')
        axs[0].legend(loc='best')
    
        axs[1].set_xlabel('k')
        axs[1].set_ylabel(f"{greek_letterz[18]}(k)")
        axs[1].legend(loc='best')
        #Add GN as gridline
            
        fig.tight_layout()
        plt.savefig(f'{top_res_path}/Figures/k_cutoff_{k_min}.svg', dpi=300, format='svg')

 

	