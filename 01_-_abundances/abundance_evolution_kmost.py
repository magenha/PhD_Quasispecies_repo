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
    
df = {}
#Upload data to DataFrame
df[43] = pd.DataFrame([])
df[30] = pd.DataFrame([])
df[43] = pd.DataFrame.from_dict(data_dict[43])
df[30] = pd.DataFrame.from_dict(data_dict[30])
df[43] = df[43].fillna(0)
df[30] = df[30].fillna(0)

temps = [30,43] 

for t in temps:

    df_t = df[t]
    
    df_t = df_t.reindex(sorted(df_t.columns), axis=1)
    L = list(df_t.columns)
    L.sort()

    k_list = [2,0,1,12]

    k_letter = {0:'B', 1:'C', 2:'A', 12:'X'}
    c_k = {2:'blue', 1:'yellow', 0:'green', 12:'black'}
    
    plt.figure()
    for i in k_list:
    
        for col in L:
            df_t[col] = df_t[col]/df_t[col].sum()
        sample = df_t.loc[i]

        x = np.array(sample.index.tolist())
        y = np.array(sample.tolist())
        

        plt.plot(x,y, 'o-', label=f"{k_letter[i]}", c=c_k[i])
    plt.ylim(bottom=3*10**-5)
    plt.yscale('log')
    plt.xlabel('step')
    plt.ylabel('relative abundance')
    plt.legend(loc='best')
    plt.title(f'Evolution of some genomes for t={t}')
    plt.savefig(f'./pics_results/abundance_evolution_kmost_{t}.svg', dpi=300, format='svg')
