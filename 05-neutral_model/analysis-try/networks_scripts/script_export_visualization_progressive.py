

import os;
import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
import pickle;
import networkx as nx;


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

top_res_path = './Topology_results/'
#top_res_path = './Topology_results/Progressive_GN/'
temps = [30, 43]
cmap = {30:'blue', 43:'red'}
symb_dict = {2:"o" , 4: "^", 10: "s", 30: "p", 60: "H"}
greek_letterz=[chr(code) for code in range(945,970)]


Abundances_path = './Data/Abundances/';
dataPathBase = './Data/'
figPathOut = './Topology_results/Progressive_GN/Figures/Visualization/'

nLoad = 1000

print('********* Searching the most abundant nodes in each step *******************')
#Get the nLoad most abundant nodes until step s at temp t
most_abundant_nodes = {t:{} for t in temps}     #most_abundant_nodes[t][step] = [list of nodes to display]

#read all files in Abundances
allFiles = os.listdir(Abundances_path)
allFiles.sort()
allFiles = [a for a in allFiles if '.csv' in a]

#First I need to know all the steps for both temps
steps = {t:[] for t in temps}

for thisFName in allFiles:
	t, step = analyze_fname(thisFName)
	steps[t].append(step)
	steps[t].sort()

for thisFName in allFiles:
	#print(thisFName)
	t, step = analyze_fname(thisFName)
	df = pd.read_csv(f'{Abundances_path}/{thisFName}', names=['hapl', 'abundance', 'relative', 'euclidean'], header=None)
	df = df.sort_values(by='abundance', ascending=False)
	k_most = df['hapl'].tolist()[:nLoad]
	if steps[t].index(step)==0:
		most_abundant_nodes[t][step] = k_most

	else:
		
		most_abundant_nodes[t][step] = list(set(most_abundant_nodes[t][steps[t][steps[t].index(step)-1]] + k_most))

print('********* Loading positions of the Grand Network *******************')

if os.path.isfile(f'{dataPathBase}GN/coordinates.co'):
    print('Coordinates already computed. Loading .......')
    with open(f'{dataPathBase}GN/coordinates.co', 'rb') as coordinates_file:
	    nodePositions = pickle.load(coordinates_file)
else:
	#Read the GN
    filepath = os.path.join(dataPathBase, "GN/edges.csv");
    thisFName = 'GrandNetwork'
    edges = np.loadtxt(filepath, delimiter=',', dtype=int);
    #Build graph
    G = nx.Graph(); 
    G.add_edges_from(edges);
    G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    #Save coordinates
    print('Computing coordinates .......')
    nodePositions = nx.spring_layout(G)
    with open(f'{dataPathBase}GN/coordinates.co', 'wb') as coordinates_file:
        pickle.dump(nodePositions, coordinates_file)


#Read a certain progressive Network
print('Loading progressive Grand Networks')
top_res_path = './Topology_results/Progressive_GN/'
allFileNames = os.listdir(f'{dataPathBase}progressive_GN/')
allFileNames.sort()

for thisFile in allFileNames:
    print(f'Processing {thisFile}')
    t,step = thisFile.split('_')
    t = int(t)
    step = int(step)
    nodes = most_abundant_nodes[t][step]
    filepath = os.path.join(dataPathBase, f"progressive_GN/{thisFile}");
    with open(filepath, 'rb') as file:
	    
        edges = pickle.load(file)
    G = nx.Graph();
    G.add_edges_from(edges);
    G = G.subgraph(nodes).copy()
    G = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    # Producing plot: 
    fig = plt.figure(); 
    #nodePositions = pos; 
    nx.drawing.nx_pylab.draw_networkx_edges(G, pos=nodePositions); 
    nx.drawing.nx_pylab.draw_networkx_nodes(G, pos=nodePositions, node_size=1);
    #print(G.nodes())
    # nx.draw(G, pos=nodePositions, node_size=abundancesToPlot); 
    fig.savefig(figPathOut +thisFile+".pdf"); 
    plt.close(fig);
    #print(filepath)
	

#plot with fixed coords

#Save fig