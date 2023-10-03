import networkx as nx
import math


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import itertools
import pandas as pd
import math
from matplotlib import animation
import random
import scipy.special
from os import listdir
from os.path import isfile, join
import os

#Read data

dir_path = '/home/samuel/Documents/PhD/Quasispecies/Data/'
index_name =  dir_path + 'seqs_index.dict'

index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])
        


WD = os.getcwd()
mypath = dir_path + 'Sequences_filtered/'

#Obtain the file names
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
data_dict = {}
data_dict[43] = {}
data_dict[30] = {}


#a = onlyfiles[0]
#Interpret file name, and extract step number from it

for a in onlyfiles:    
    if a[0] == 'c':
        step = int(a.split('-')[1][1:])
        t=43
        print(f'step={step} for t={t}')
    elif a[0] == '3':
        step = int(a.split('-')[2])
        t=30
        print(f'step={step} for t={t}')
    
    file_name = dir_path + 'Sequences_filtered/' + a
    seq_2_ab = {}
    with open(file_name, 'r') as r:
        for line in r:
            if line[0] == '>':
                abundance = int(line.split('-')[1][:-1])
                #print(abundance)
            else:
                sequence = line[:-1]
                hapl = index_dict[sequence]
                seq_2_ab[hapl]= abundance
    data_dict[t][step] = seq_2_ab


seqs_list = list(index_dict.keys())

b = seqs_list[0]
c = seqs_list[1]
a = seqs_list[2]
x = seqs_list[12]

d_ab=0
d_ac=0
d_bc=0
d_ax = 0
d_bx = 0
d_cx = 0
for char in range(len(a)):
    if a[char] != b[char]:
        print(f'Diference a -> b: {char+1} {a[char]} {b[char]}')
        d_ab+=1
    if a[char] != c[char]:
        print(f'Diference a -> c: {char+1} {a[char]} {c[char]}')
        d_ac+=1
    if b[char] != c[char]:
        print(f'Diference c -> b: {char+1} {c[char]} {b[char]}')
        d_bc+=1
    if b[char] != x[char]:
        print(f'Diference b -> x: {char+1} {b[char]} {x[char]}')
        d_bx+=1
    if a[char] != x[char]:
        print(f'Diference a -> x: {char+1} {a[char]} {x[char]}')
        d_ax+=1



#Plotting

a=1
G=nx.Graph()
G.add_edges_from([(1,2),(2,3),(3,4),(4,1),
                  (5,1), (5,2), (6,2), (6,3),
                 (7,3), (7,4),(8,1), (8,4),
                 (9,1),(10,2),(11,3), (12,4)]) #define G
fixed_positions = {1:(-a/2,a/2),
                   2:(-a/2,-a/2),
                  3: (a/2,-a/2),
                  4: (a/2,a/2),
                  5: (-a,0),
                  6: (0,-a),
                  7: (a,0),
                  8: (0,a),
                  9: (-a*math.sqrt(2),a*math.sqrt(2)),
                  10: (-a*math.sqrt(2),-a*math.sqrt(2)),
                  11: (a*math.sqrt(2),-a*math.sqrt(2)),
                  12: (a*math.sqrt(2),a*math.sqrt(2))}#dict with two of the positions set
nodes_names = {1:"A", 2:"B", 3:"C", 4:"X", 5:"A&B", 6:"B&C", 7:"C&X", 8:"A&X", 
               9:"A others", 10:"B others", 11:"C others", 12:"X others"}
fixed_nodes = fixed_positions.keys()
pos = nx.spring_layout(G,pos=fixed_positions, fixed = fixed_nodes)
#G = nx.relabel_nodes(G, nodes_names)
nx.draw_networkx(G,pos, labels=nodes_names, with_labels=True)
#nx.draw_networkx_labels(G,pos, labels=nodes_names,ax=ax[0], with_labels=True)
#ax[1].set_ylim(tuple(i*1.1 for i in ax[1].get_ylim()))
plt.show()