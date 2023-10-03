'''
    This script is intended to check if we have been generating properly the progressive GN
'''

import pickle;
import networkx as nx;
import os;
import numpy as np;
import matplotlib.pyplot as plt;

#Sum of the final GNs at each temp == GN (Luis)

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

#Upload all edges. 

datafile = './Data/progressive_GN/'

temps = [30,43]
cmap = {30:'blue', 43:'red'}

#Activate this to collapse both temps GNs
if 0==1:
    file_30 = datafile + '30_60'
    file_43 = datafile + '43_60'




    edges = {t:[] for t in temps}

    for t in temps:
        with open(datafile + f'{t}_60', 'rb') as file:
            edges[t] = pickle.load(file)

    print('Loaded GN at each temp')
    edg_list = edges[30]
    edges2 = edges[43]
    print('Mixing')
    i=0
    for edge2 in edges2:
        #print(edge2, edge2[0], edge2[1])
        if edge2 in edg_list: pass
        elif edge2[::-1] in edg_list: pass
        else: edg_list.append(edge2)
        if i%10000 == 0:
            print(f'{round(i/len(edges2)*100, 2)}%')
        i+=1
    with open('./Data/progressive_GN/test/sumGN', 'wb') as file:
        pickle.dump(edg_list, file)

    print([len(edges[t]) for t in temps], len(edg_list))



#Activate this to know the properties of each final GN at each temp
if 0==1:
    file_30 = datafile + '30_60'
    file_43 = datafile + '43_60'

    G = {t: nx.Graph() for t in temps}
    edges = {t:[] for t in temps}

    for t in temps:
        with open(datafile + f'{t}_60', 'rb') as file:
            edges[t] = pickle.load(file)
    
            G[t].add_edges_from(edges[t])
            print(f'For t={t} we have', len(G[t].nodes()), 'nodes and ', len(G[t].edges()), 'edges')



#Activate this to check if the nNodes in each subnetwork is the same as the number of haplotypes in each step
if 0==1:
    dataPathBase = '/home/samuel/Documents/samuel/Quasispecies/04_-_networks/work_withoutN/Reconstruction_QBeta/Data/'
    # Reading number of haplotypes: 
    fIn = open(os.path.join(dataPathBase+"fileNames.csv"), 'r'); 
    allFileNames = fIn.read().splitlines(); 
    fIn.close(); 
    allFileNames.sort()

    haplotypes = {t:{} for t in temps}
    nodes = {t:{} for t in temps}

    for thisFile in allFileNames:
        i=0
        t, step = analyze_fname(thisFile)
        filepath = os.path.join(dataPathBase, thisFile + '.fasta');
        with open(filepath, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line[0] == '>':
                    pass
                else:
                    i+=1
        haplotypes[t][step] = i


    #Reading number of nodes in each subnetwork
    subNetworkPath = os.path.join(dataPathBase, "GN/SubNetworks")
    for thisFName in allFileNames:
        t,step = analyze_fname(thisFName);
        print(f'Analyzing {thisFName}')
        #Loading edges
        filepath = os.path.join(subNetworkPath, "edges_" + thisFName +".csv");
        edges = np.loadtxt(filepath, delimiter=',', dtype=int); 
        #Build graph
        G = nx.Graph(); 
        G.add_edges_from(edges);

        nodes[t][step] = len(list(G.nodes()))
            
    for t in temps:
        L = list(haplotypes[t].keys())
        L.sort()
        x = [nodes[t][step] for step in L]
        y = [haplotypes[t][step] for step in L]
        plt.scatter(x,y, c=cmap[t])
        #for step in L:
        #    print(t, step, haplotypes[t][step], nodes[t][step])
    t=43
    sample = haplotypes[t].values()
    sample = np.linspace(min(sample)*0.8, max(sample)*1.1, 100)
    plt.plot(sample, sample, '--', c='black')
    plt.xlabel('nodes in subnetwork')
    plt.ylabel('number of haplotypes')
    plt.savefig('/home/samuel/Documents/samuel/Quasispecies/04_-_networks/work_withoutN/Reconstruction_QBeta/Topology_results/nodes-haplotypes.svg', dpi=300, format='svg')
    plt.show()


















#Activate this to compare directly the sum of progGN with final GN
if 0==1:
    #Read the added GN
    with open('./Data/progressive_GN/test/sumGN', 'rb') as file:
        GN_added_edges = pickle.load(file)
    #Read LFSeoane's GN

    print('Analyzing final Grand Network')
    dataPathBase = "./Data/"; 
    filepath = os.path.join(dataPathBase, "GN/edges.csv");
    thisFName = 'GrandNetwork'
    GN_edges = np.loadtxt(filepath, delimiter=',', dtype=int);
    #Build graph
    GN = nx.Graph(); 
    GN.add_edges_from(GN_edges);

    G_added = nx.Graph()
    G_added.add_edges_from(GN_added_edges)

    #Check if all edges in G_added are in GN
    print('Checking if all edges in summed_GN are in computed GN')
    ok=0
    _ok = 0
    for edge in GN_added_edges[:10000]:
        if edge in GN_edges:
            ok+=1
        else:
            _ok +=1
    print('yes:', ok, 'not:', _ok)