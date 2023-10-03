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
cmap = {30:'b', 43:'r'}
out_dir = '/home/samuel/Documents/PhD/Quasispecies/Quasispecies_evolution/01_-_abundances'
n_max = 15
r_max = 10000
# Get the number of sequences obtained in each step

plt.figure()
for t in temps:
    datadict = {}
    L = df[t].columns.tolist()
    L.sort()
    for element in L:
        datadict[element] = int(df[t][element].sum())
    
    dataf = pd.DataFrame.from_dict(datadict, orient='index', columns=['Total'])
    #dataf.sort_values('Total', ascending=True)
    plt.plot(dataf.index.tolist(), dataf['Total'].tolist(), 'o-', label=f't={t}', c=cmap[t])


plt.ylabel('total of sequences')
plt.xlabel('step')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(f'{out_dir }/pics_results/total_sequences_Nremrem.svg', dpi=300, format='svg')

plt.figure()
#Haplotypes
for t in temps:
    datadict = {}
    L = df[t].columns.tolist()
    L.sort()
    for element in L:
        df_x = df[t][element]
        datadict[element] = int(len(df[t][element].to_numpy().nonzero()[0]))
    
    dataf = pd.DataFrame.from_dict(datadict, orient='index', columns=['Total'])
    #dataf.sort_values('Total', ascending=True)
    plt.plot(dataf.index.tolist(), dataf['Total'].tolist(), 'o-', label=f't={t}', c=cmap[t])

    
plt.ylabel('total of haplotypes')
plt.xlabel('step')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(f'{out_dir }/pics_results/total_haplotypes_Nremrem.svg', dpi=300, format='svg')

plt.figure()
# number of haplotypes vs. number of sequences

for t in temps:
    datadict_htyp = {}
    datadict_nseq = {}
    L = df[t].columns.tolist()
    L.sort()
    for element in L:
        #df_x = df[t][element]
        datadict_htyp[element] = int(len(df[t][element].to_numpy().nonzero()[0]))
        datadict_nseq[element] = int(df[t][element].sum())
    
    dataf_htyp = pd.DataFrame.from_dict(datadict_htyp, orient='index', columns=['Total'])
    dataf_nseq = pd.DataFrame.from_dict(datadict_nseq, orient='index', columns=['Total'])
    #dataf.sort_values('Total', ascending=True)
    xs = dataf_nseq['Total'].tolist()
    ys = dataf_htyp['Total'].tolist()
    plt.plot(xs, ys, 'o', label=f't={t}', c=cmap[t])
    i=0
    for x,y in zip(xs,ys):

        label = "{:d}".format(L[i])

        plt.annotate(label, # this is the text
                    (x,y), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,10), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
        i +=1
plt.ylabel('haplotypes')
plt.xlabel('sequences')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(f'{out_dir }/pics_results/seqs_hapls_Nremrem.svg', dpi=300, format='svg')