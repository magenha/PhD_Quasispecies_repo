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

dir_path = './Data-lethal/'
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
n_max = 15
r_max = 10000

for t in temps:

    df_t = df[t]

    L = list(df_t.columns)
    L.sort()
    n = len(L)
    

    abundances = {i: [] for i in L}

    for i,step in enumerate(L):        

        vlist = df_t[step].to_list()
        vlist.sort(reverse=True)
        vector = np.array(vlist)
        # Normalize by 1
        vector = vector / vector.sum()

        abundances[step]= vector
    
    
    
    for i,step in enumerate(L):
        plt.figure()
        x = [1+i for i in range(len(abundances[step]))]
        plt.plot(x, abundances[step], 'o-', label=f'{L[i]}')
            

        plt.yscale('log')
        plt.xscale('log')
        plt.ylabel('abundance')
        plt.xlabel('r')
        plt.ylim(ymax=0.8, ymin=5*10**-7)
        plt.xlim(xmax=2*10**6, xmin=0.5)
        plt.title(f't={t}')
        plt.legend()
        plt.savefig(f'./temp/sample_{i}_{t}.png')

for t in temps:
    # Create new figure for GIF
    fig, ax = plt.subplots()

    # Adjust figure so GIF does not have extra whitespace
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    ax.axis('off')
    ims = []

    for m in range(n):
        im = ax.imshow(plt.imread(f'./temp/sample_{m}_{t}.png'), animated = True)
        ims.append([im])

    speed = 400
    ani = animation.ArtistAnimation(fig, ims, interval=speed)
    ani.save(f'./pics_results/evol_abundances_{t}.gif')