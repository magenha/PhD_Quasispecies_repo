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

dir_path = './Data/'
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
n_max = 10
r_max = 10000

for t in temps:

    df_t = df[t]
    
    #Set a fixed rank
    df_t = df_t.sort_values(by=2, ascending=False)
    
    
    L = list(df_t.columns)
    L.sort()
    
    r_max = 2**n_max

    x = [2**i for i in range(n_max+1)]
    z = {i: {j: 0 for j in range(L[-1]+1)} for i in range(n_max+1)}

    for k,step in enumerate(L):


        vector = df_t[step]
        vector_index = [i for i in range(r_max)]
        vector = vector / vector.sum()
        vector_list = vector.tolist()
        

        for j in vector_index:
            if j ==0:
                i_sampl = 0
                z[i_sampl][step] += vector_list[j] 
            else: 
                i_sampl = int(math.log2(j))
                print(i_sampl, step)
                print(vector_list[j])
                z[i_sampl][step] += vector_list[j] 
        
        for i in range(n_max+1):
            z[i][step] = np.log10(1.0+z[i][step])
        
    for k,step in enumerate(L):
        if k<len(L)-1:
            d_s = L[k+1]-L[k]
            if d_s !=1:
                for m in range(1,d_s):
                    for key in z.keys():
                        z[key][L[k]+m] = 0.5*z[key][L[k]] + 0.5*z[key][L[k+1]]

        else:
            pass



    img = np.zeros((L[-1]+1,n_max+1))

    for step in range(max(L)+1):
        img[step] = np.array([z[x][step] for x in range(n_max+1)])

    #Fill empty initial steps
    end_zero = 0
    for i in range(max(L)+1):
        
        if img[i].sum() != 0:
            break
        else:
            end_zero +=1
    for i in reversed(range(end_zero)):
        img[i] = img[end_zero]
    img = img.T
    

    plt.figure()
    plt.pcolormesh(img)
    plt.set_cmap('hot')
    #plt.xscale('symlog')
    plt.xlabel('   ')
    plt.ylabel('   ')
    
    #y_labels = [L[i] for i in range(len(L))]
    #plt.yticks(y_labels)
    plt.title(f't={t}')
    plt.tight_layout()
    
    #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    cbar = plt.colorbar(location='right')
    #cbar.ax.set_ylabel('log(Abundance)', rotation=270)
    plt.savefig(f'./abundances_evolution_heatmap_fixrank_{t}.svg', dpi=300, format='svg')
'''
for t in temps:

    df_t = df[t]
    
    #Set a fixed rank
    df_t = df_t.sort_values(by=2, ascending=False)
    
    
    L = list(df_t.columns)
    L.sort()
    r_max = 2**n_max

    x = [2**i for i in range(n_max+1)]
    z = {i: {j: 0 for j in range(L[-1]+1)} for i in range(n_max+1)}

    for k,step in enumerate(L):


        vector = df_t[step]
        vector_index = [i for i in range(r_max)]
        vector = vector / vector.sum()
        
        vector = vector.sort_values(ascending=False)
        vector_list = vector.tolist()

        #z[j] = #Abundances in step=2 binned with x


        for j in vector_index:
            if j ==0:
                i_sampl = 0
                z[i_sampl][step] += vector_list[j] 
            else: 
                i_sampl = int(math.log2(j))

                z[i_sampl][step] += vector_list[j] 
            
        for i in range(n_max+1):
            z[i][step] = np.log10(1.0+z[i][step])
        
    for k,step in enumerate(L):
        if k<len(L)-1:
            d_s = L[k+1]-L[k]
            if d_s !=1:
                for m in range(1,d_s):
                    for key in z.keys():
                        #z[key][L[k]+m] = 1.0*z[key][L[k]] + 0.0*z[key][L[k+1]]
                        z[key][L[k]+m] = 0.5*z[key][L[k]] + 0.5*z[key][L[k+1]]
                        #z[key][L[k]+m+1] = ((1+m)/d_s)*(z[key][L[k]])+(d_s-(1+m))/d_s*(z[key][L[k+1]])

        else:
            pass
        
    #for i in range(n_max+1):
    #    for step in range(L[-1]):
    #        z[i][step] = np.log2(1.0+z[i][step])


    img = np.zeros((L[-1]+1,n_max+1))
    #img = np.zeros((len(L),n_max+1))
    for step in range(max(L)+1):
        img[step] = np.array([z[x][step] for x in range(n_max+1)])
        #-1.0*np.log10(1.0+np.array(abundances[i]))
    #Fill empty steps
    end_zero = 0
    for i in range(max(L)+1):
        
        if img[i].sum() != 0:
            break
        else:
            end_zero +=1
    for i in reversed(range(end_zero)):
        img[i] = img[end_zero]
        
    
    #Improve the x axis

    plt.figure()
    plt.pcolormesh(img)
    plt.set_cmap('hot')
    #plt.xscale('symlog')
    plt.ylabel('step')
    plt.xlabel('Rank (log2)')
    
    #y_labels = [L[i] for i in range(len(L))]
    #plt.yticks(y_labels)
    plt.title(f't={t}')
    plt.tight_layout()
    
    #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    #, pad=0.15
    plt.colorbar(location='right')
    plt.savefig(f'./pics_results/abundances_evolution_heatmap_{t}')'''