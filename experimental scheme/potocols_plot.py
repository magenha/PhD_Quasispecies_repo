#import libraries

from os import listdir
from os.path import isfile, join
import os
import itertools
import random
import time
from datetime import datetime
random.seed(datetime.now().timestamp())


import numpy as np
import pandas as pd
import math


import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d



#Read and storage data

#Read the index dictionary
#Revert the index 

print('Reading seqs_index', datetime.now())
data_path = '/home/samuel/Documents/PhD/Quasispecies/Data'
index_name =  data_path+'/seqs_index.dict'

index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])
        
print('Loading data', datetime.now())

mypath = data_path + '/Sequences_filtered/'


#Obtain the file names
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
data_dict = {}
data_dict[43] = {}
data_dict[30] = {}

#Interpret file name, and extract step number from it

for a in onlyfiles:    
    if a[0] == 'c':
        step = int(a.split('-')[1][1:])
        t=43
        #print(f'step={step} for t={t}')
    elif a[0] == '3':
        step = int(a.split('-')[2])
        t=30
        #print(f'step={step} for t={t}')
    
    file_name = mypath + a
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

df = {}
#Upload data to DataFrame
df[43] = pd.DataFrame([])
df[30] = pd.DataFrame([])
df[43] = pd.DataFrame.from_dict(data_dict[43])
df[30] = pd.DataFrame.from_dict(data_dict[30])
df[43] = df[43].fillna(0)
df[30] = df[30].fillna(0)

temps = [30,43] 
cmap = {30:'b', 43:'r'}

plt.figure()
for t in temps:
    df_t = df[t]
    L=list(df_t.columns)
    L.sort()
    y = [t for i in range(len(L))]
    #x = [1] + L
    #y = [37] + [t for i in range(len(L)-1)]
    #plt.plot([1],[37], c=cmap[t])
    L = list(itertools.chain([0], L))
    y = list(itertools.chain([37], y))
    plt.plot(L,y, 'o-', c=cmap[t])
out_dir = '/home/samuel/Documents/PhD/Quasispecies/Quasispecies_evolution/experimental scheme/'
plt.ylabel(f't')
plt.xlabel(f'steps')
plt.savefig(out_dir+'protocol.svg', dpi=300, format='svg')
