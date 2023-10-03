
'''
    This script is intended to study the difussion dynamics of the Qbeta virus. 

'''

import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

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

def calculateParametersLSR(x, y):
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



	return b, w, SEw


#Select the wt
wt_sequence='CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG'
temps = [30,43]
cmap = {30: 'blue', 43:'red'}
nLoad = 1000

def distance(seq, seq_ref):
    '''
        Function that returns the Humming distance between
        two sequences
        Input:
            -sequence to compare
            -reference sequence
        Return:
            -int : distance
    '''
    d=0
    for char in range(len(seq_ref)):
        if seq[char] != seq_ref[char]:
            d+=1
    return d

dir_path = '/home/samuel/Documents/PhD/Quasispecies/Data/'
index_name =  dir_path + 'seqs_index.dict'

index_dict = {}
inv_index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])
        inv_index_dict[int(L[0])] = L[1][:-1]

#print(list(inv_index_dict.keys())[:10])

datapath = '/home/samuel/Documents/PhD/Quasispecies/Data/Sequences_filtered/N_rem_rem/'
#Read all fileNames in folder
allFileNames = [f for f in listdir(datapath) if isfile(join(datapath, f))]
allFileNames.sort()

d_results_list = {t:{} for t in temps}

for thisFileName in allFileNames:
    print(f'Reading {thisFileName}')
    t,step = analyze_fname(thisFileName)
    d_results_list[t][step] = []

    with open(os.path.join(datapath, thisFileName), 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '>':
                id = int(line.split('-')[0].replace('>', ''))
                abundance=int(line.split('-')[1].replace('\n', ''))
                #print(id)
                d = distance(inv_index_dict[id],wt_sequence)
                for n in range(abundance):
                    d_results_list[t][step].append(d)
'''
t=30
step=10
print(d_results_list[t][step])
bins = np.arange(0,20, 1)
print(np.array(d_results_list[t][step]).mean())
print(np.array(d_results_list[t][step]).std())
plt.hist(d_results_list[t][step], bins=bins)
plt.show()'''

#Active this to plot mean values
if 1==1:
    for t in temps:
        L = list(d_results_list[t].keys())
        L.sort()
        y = []
        y_err = []
        y_sq = []
        index = []
        step_list=[]

        #Compute mean d
        for step in L:
            if step: #not in  [1,4,10,60]
                step_list.append(step)
                mean = np.array(d_results_list[t][step]).mean()
                std = np.array(d_results_list[t][step]).std()
                y_sq.append((np.array(d_results_list[t][step])**2).mean())
                y.append(mean)
                y_err.append(std**2)
                index.append((std**2)/mean)
        #plt.plot(step_list, y_err,'o-', c=cmap[t])
        #plt.ylabel('variance')
        #plt.xlabel('step')
        #plt.tight_layout()
        #plt.savefig(f'./difussion_variance.svg', dpi=300, format='svg')
        #plt.figure()
        b,w,E = calculateParametersLSR(np.log2(np.array(step_list)), np.log2(np.array(y_sq)))
        plt.plot(step_list, y_sq,  'o-', c=cmap[t], label=f'{t}')

        x = np.linspace(0.9, 61, 100)
        y_aster = 2**b *x**w
        plt.plot(x,y_aster, '--', c='black', label=f'm={round(w,3)}+-{round(E,3)}')

    plt.legend()
    plt.ylabel('<r^2>')
    plt.xlabel('step')
    plt.tight_layout()
    #plt.yscale('log')
    #plt.xscale('log')
    plt.savefig(f'./MSD_fit.svg', dpi=300, format='svg')
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(f'./MSD_fit_log.svg', dpi=300, format='svg')
    plt.show()


#Read each temp and step

#Find the distribution of d_H respect to wt

#save mean and std dev

#Plot mean and std vs time fot both temps

