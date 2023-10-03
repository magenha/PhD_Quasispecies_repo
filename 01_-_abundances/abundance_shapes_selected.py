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
mypath = dir_path + 'Sequences_filtered/N_rem_rem/'

#Obtain the file names
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#onlyfiles = [a for a in onlyfiles if '_Nrem' in a]
#onlyfiles = [a for a in onlyfiles if '_Nrem' not in a]
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
    
    file_name = dir_path + 'Sequences_filtered/N_rem_rem/' + a
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
n_max = 15
r_max = 10000

#paramters for fit
m = -1.22
c = -2.0

for t in temps:

    df_t = df[t]

    L = list(df_t.columns)
    L.sort()
    n = len(L)
    

    abundances = {i: [] for i in L}

    m_list = []
    c_list = []
    r_list = []

    for i,step in enumerate(L):        

        vlist = df_t[step].to_list()[:r_max]
        vlist.sort(reverse=True)
        vector = np.array(vlist)
        # Normalize by 1
        vector = vector / vector.sum()

        abundances[step]= vector
    
    x = [1+i for i in range(r_max)]
    plt.figure()
    for i,step in enumerate(L):
        if step in [8, 10, 30, 40]:
            plt.plot(x, abundances[step], 'o-', label=f'{L[i]}')

    '''
    for element in L:
        x = []
        y = []
        for i,el in enumerate(abundances[element]):
            if el == 0.0:
                pass
            else:
                x.append(i+1)
                y.append(el)

        y_log = np.log(y)
        x_log = np.log(np.array(x))
        A = np.vstack([x_log, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y_log, rcond=None)[0]
        m_list.append(m)
        c_list.append(c)
        r = np.corrcoef(x_log, y_log)
        
        r_list.append(r[0,1]**2)
    m = np.mean(np.array(m_list))
    c = np.mean(np.array(c_list))
    r_sq = np.mean(np.array(r_list))
    
    #Plot fitted power-law
    x = np.arange(0,r_max+1, dtype=int)
    y = np.e**(c)*np.power(x,m)
    plt.plot(x,y,'--', c='brown', label='fitting', linewidth=5)  
    print(f't={t} fitting m={m}, c={c}, R^2={r_sq}') 
    '''
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('abundance')
    plt.xlabel('r')
    plt.ylim(ymax=0.8, ymin=5*10**-7)
    plt.title(f't={t}')
    plt.legend()
    plt.savefig(f'./pics_results/svg/abundances_shape_{t}_without_Nremrem.svg', dpi=300, format='svg')

