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
import pickle

#Directory holding all data
dir_path = '/home/samuel/Documents/PhD/Quasispecies/Quasispecies_evolution/05-neutral_model/analysis-try/Data/'

#Name of the index file
index_name =  dir_path + 'seqs_index.dict'
#Folder contaning the files
mypath = dir_path 

#Output directory
out_dir = '/home/samuel/Documents/PhD/Quasispecies/Quasispecies_evolution/05-neutral_model/analysis-try/overlapping/Results/'

#############################################################################

#Read index
index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])
        

#Read files if this is the first time

file_to_save = out_dir + 'Computations/df.pickle'
try:
    
    with open(file_to_save, 'rb') as f:
        df = pickle.load(f)
    print('df already saved. Reading...')
except:
    print('Computing df ...')

    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyfiles = [a for a in onlyfiles if 'all_trim_merged_filter' in a]
    print(onlyfiles)

    data_dict = {}
    data_dict[43] = {}
    data_dict[30] = {}

    for a in onlyfiles:    
        if a[0] == 'c':
            step = int(a.split('-')[1][1:])
            t=43
            print(f'step={step} for t={t}')
        elif a[0] == '3':
            step = int(a.split('-')[2])
            t=30
            print(f'step={step} for t={t}')
        
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

    with open(file_to_save, 'wb') as f:
        pickle.dump(df,f)


temps = [30,43] 
cmap = {30:'b', 43:'r'}

#Compute phi and theta
measures = {m:{t:{} for t in temps} for m in ['phi', 'theta']}

for t in temps:
    
    try:
        with open(f'{out_dir}Computations/phi_{t}', 'rb') as f:
            measures['phi'][t]['normal'] = pickle.load(f)
        with open(f'{out_dir}Computations/phi_{t}_log', 'rb') as f:
            measures['phi'][t]['log'] = pickle.load(f)
        with open(f'{out_dir}Computations/th_{t}', 'rb') as f:
            measures['theta'][t]['normal'] = pickle.load(f)
        with open(f'{out_dir}Computations/th_{t}_log', 'rb') as f:
            measures['theta'][t]['log'] = pickle.load(f)
    
        
        print(f'Loaded phi and theta for t={t}')
    except:
        print(f'Computing phi and theta for t={t}')

        df_t = df[t]
        df_log = (df_t + 1).apply('log')
        L = df_t.columns
        L = list(L)
        L.sort()
        
        df_norm = pd.DataFrame([])
        df_norm.index=df_t.index
        df_norm_log = pd.DataFrame([])
        df_norm_log.index=df_t.index

        for x in L:
            arr = np.array(df_t[x].tolist())
            normalized_arr = arr / np.sqrt(np.sum(arr**2))
            df_norm[x] = normalized_arr
            
            arr = np.array(df_log[x].tolist())
            normalized_arr = arr / np.sqrt(np.sum(arr**2))
            df_norm_log[x] = normalized_arr

        cos_phi = []
        for i in range(len(L)-1):
            a = np.array(df_norm[L[i]].tolist()) 
            b = np.array(df_norm[L[i+1]].tolist()) 
            cos_phi.append(np.dot(a,b))
        cos_phi = np.around(cos_phi, decimals=6)
        phi = np.arccos(np.array(cos_phi))*180/np.pi
        
        with open(f'{out_dir}Computations/phi_{t}', 'wb') as f:
            pickle.dump(phi,f)

        cos_phi_log = []
        for i in range(len(L)-1):
            a = np.array(df_norm_log[L[i]].tolist()) 
            b = np.array(df_norm_log[L[i+1]].tolist()) 
            cos_phi_log.append(np.dot(a,b))
        cos_phi_log = np.around(cos_phi_log, decimals=6)
        phi_log = np.arccos(np.array(cos_phi_log))*180/np.pi
        with open(f'{out_dir}Computations/phi_{t}_log', 'wb') as f:
            pickle.dump(phi_log,f)

        #Compute theta from the origin step
        #theta

        cos_th = []
        cos_th_log = []
        origin = np.array(df_norm[L[0]].tolist())
        origin_log = np.array(df_norm_log[L[0]].tolist())

        for i in range(len(L)-1):
            cos_th.append(np.dot(origin, np.array(df_norm[L[i]].tolist())))
        cos_th = np.around(cos_th, decimals=6)
        th = np.arccos(np.array(cos_th))*180/np.pi
        
        with open(f'{out_dir}Computations/th_{t}', 'wb') as f:
            pickle.dump(th,f)

        for i in range(len(L)-1):
            cos_th_log.append(np.dot(origin_log, np.array(df_norm_log[L[i]].tolist())))
        cos_th_log = np.around(cos_th_log, decimals=6)
        th_log = np.arccos(np.array(cos_th_log))*180/np.pi
        
        with open(f'{out_dir}Computations/th_{t}_log', 'wb') as f:
            pickle.dump(th_log,f)
        
        measures['phi'][t]['normal'] = phi
        measures['phi'][t]['log'] = phi_log
        measures['theta'][t]['normal'] = th
        measures['theta'][t]['log'] = th_log




#plots

#plot phi
plt.figure()
for t in temps:
    
    x = list(df[t].columns)
    x.sort()
    y = [measures['phi'][t]['normal'][i]/(x[i+1]-x[i]) for i in range(0,len(x)-1)]
    plt.plot(x[:-1], y ,'o-' , c=cmap[t],label=f'{t}')
plt.xlabel('step')
plt.ylabel('phi')
plt.legend(loc='best')
plt.savefig(f'{out_dir}Images/phi_normal.svg', dpi=300, format='svg')


#plot theta
plt.figure()
for t in temps:
    x = list(df[t].columns)
    x.sort()
    plt.plot(x[1:],measures['theta'][t]['normal'], 'o-' , c=cmap[t],label=f'{t}')
plt.xlabel('step')
plt.ylabel('theta')
plt.legend(loc='best')
plt.savefig(f'{out_dir}Images/theta_normal.svg', dpi=300, format='svg')

#plot phi log
plt.figure()
for t in temps:
    #print(f'philog {t}', measures['phi'][t]['log'])
    x = list(df[t].columns)
    x.sort()
    y = [measures['phi'][t]['log'][i]/(x[i+1]-x[i]) for i in range(0,len(x)-1)]
    plt.plot(x[:-1], y, 'o-' , c=cmap[t],label=f'{t}')
    #plt.plot(x[:-1],measures['phi'][t]['log'], 'o-' , c=cmap[t],label=f'{t}')
plt.xlabel('step')
plt.ylabel('phi log')
plt.legend(loc='best')
plt.savefig(f'{out_dir}Images/phi_log.svg', dpi=300, format='svg')


#plot theta log
plt.figure()
for t in temps:
    x = list(df[t].columns)
    x.sort()
    plt.plot(x[1:],measures['theta'][t]['log'], 'o-' , c=cmap[t],label=f'{t}')
plt.xlabel('step')
plt.ylabel('theta log')
plt.legend(loc='best')
plt.savefig(f'{out_dir}Images/theta_log.svg', dpi=300, format='svg')

#Plot sumphi
aspects = ['normal', 'log']
for aspect in aspects:
    data = {t:{} for t in temps}
    for t in temps:
        x = list(df[t].columns)
        x.sort()
        phi = measures['phi'][t][aspect]
        #print(t,phi)
        y_phi = [phi[0]]
        for i in range(1,len(phi)):
            y_phi.append(y_phi[-1] + phi[i])
            
        y_phi = np.array(y_phi)
        #print(t, y_phi)
        with open(f'{out_dir}Computations/cumphi_{t}_{aspect}', 'w') as f:
            for i in range(len(y_phi)):
                f.write(f'{x[i]}, {y_phi[i]}\n')
        with open(f'{out_dir}Computations/cumphi_{t}_{aspect}', 'r') as f:
            lines = f.readlines()
            for line in lines:
                x,y = line[:-1].split(', ')
                data[t][int(x)] = float(y)
        
    plt.figure()
    for t in temps:
        x=[]
        y=[]
        x = list(data[t].keys())
        y = np.array([data[t][step] for step in data[t].keys()])
        
        plt.plot(x,y, 'o-',c=cmap[t] ,label=f't={t}')


    plt.xlabel('step')
    plt.ylabel('sum(phi)')
    plt.legend(loc='best')
    plt.yscale('linear')
    plt.xscale('linear')
    plt.savefig(f'{out_dir}Images/cumphi_{aspect}.svg', dpi=300, format='svg')
