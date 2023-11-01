
'''
    This script is intended to study the difussion dynamics of the Qbeta virus. 

'''

import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FuncFormatter

def format_func(value, tick_number):
    return f'{value:.2f}'

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

greek_leeterz = [chr(code) for code in range(945,970)]
#Select the wt
wt_sequence_A = 'CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG'
wt_sequence_B='CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCATAAAGCGTTGAAACTTTG'
wt_sequence_C = 'CAACAAGGTCAGCTATATCATAATATCGGTATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCATAAAGCGTTGAAACTTTG'
wt_sequence_D = 'CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGACTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG'

wt_sequence = wt_sequence_D

temps = [30,43]
cmap = {30: '#42c1eb', 43:'#ed4242'}
cmap_line = {30: 'blue', 43:'red'}
symbmap = {30:'o', 43:'^'}
symbdict = {2:'o',12:'s', 30:'^', 58:'*'}
symbline = {30: '--', 43: '-'}
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
        for i,line in enumerate(lines):
            if line[0] == '>':
                id = int(line.split('-')[0].replace('>', ''))
                abundance=int(line.split('-')[1].replace('\n', ''))
                sequence = lines[i+1].replace('\n', '')
                #print(id)
                d = distance(sequence,wt_sequence)
                for n in range(abundance):
                    #d_results_list[t][step].append(math.sqrt(d))
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

#Active this to plot d^2 dependence with time
if 1==1:
    plt.figure(figsize=(8,6))
    for t in temps:
        #for t in [30]:
        L = list(d_results_list[t].keys())
        L.sort()
        y = []
        y_err = []
        y_sq = []
        index = []
        step_list=[]

        #Compute moments of d
        for step in L:
            if step not in  [1,4,10,60] and step >30:
                step_list.append(step)

                y_sq.append((np.array(d_results_list[t][step])**2).mean())

        b,w,E = calculateParametersLSR(np.log2(np.array(step_list)), np.log2(np.array(y_sq)))
        
        
        x = np.array(step_list)
        y_aster = 2**b *x**w
        plt.plot(x,y_aster, symbline[t], c=cmap_line[t], label=f'{greek_leeterz[0]}:{round(w,2)}+-{round(E,2)}')
        step_list = []
        y_sq = []

        for step in L:
            if step not in  [1,4,10,60] and step <=30:
                step_list.append(step)

                y_sq.append((np.array(d_results_list[t][step])**2).mean())

        #Do fitting for steps <25
        #step_list = np.array(step_list)
        #b,w,E = calculateParametersLSR(np.log2(np.array(step_list)), np.log2(np.array(y_sq)))
        
        
        #x = np.array(step_list)
        #y_aster = 2**b *x**w

        #plt.plot(x,y_aster, '-.', c=cmap[t], label=f'<30 {greek_leeterz[0]}:{round(w,2)}+-{round(E,2)}')

        step_list = []
        y_sq = []
        #Now search the points to do the plot of r^2
        for step in L:
            if step not in [1,4,10,60]:
                step_list.append(step)

                y_sq.append((np.array(d_results_list[t][step])**2).mean())

        plt.plot(np.array(step_list), np.array(y_sq), symbmap[t], c=cmap[t], label=f'{t}')
    
    plt.yticks(np.linspace(min(y_sq), max(y_sq)*1.1, 5))
    plt.xticks(np.linspace(min(step_list), max(step_list), 5))
    plt.rcParams.update({'font.size': 18})
    plt.gca().yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.ylabel('Average number of mutations, m')
    plt.xlabel('step')
    plt.tight_layout()
    plt.savefig(f'./MSD_fit.svg', dpi=300, format='svg')
    plt.savefig(f'./MSD_fit.png', dpi=300, format='png')
    
    #plt.yticks(np.linspace(min(y_sq), max(y_sq)*1.1, 5))
    #plt.xticks(np.linspace(min(step_list), max(step_list), 5))
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.xlim((28,61))
    #plt.rcParams.update({'font.size': 18})
    #plt.gca().yaxis.set_major_formatter(FuncFormatter(format_func))
    #plt.savefig(f'./MSD_fit_log.svg', dpi=300, format='svg')
    #plt.savefig(f'./MSD_fit_log.png', dpi=300, format='png')

#Activate this to plot <r>, variance, index with fitting
if 1==1:
    #Mean
    plt.figure(figsize=(8,6))
    for t in temps:
        #for t in [30]:
        L = list(d_results_list[t].keys())
        L.sort()
        y = []
        step_list=[]

        #Compute moments of d
        for step in L:
            if step not in [1,4,10,60]:
                step_list.append(step)
                mean = np.array(d_results_list[t][step]).mean()
                y.append(mean)
        y = np.array(y)
        #plt.plot(np.array(step_list), np.array(y), 'o-', c=cmap[t], label=f'{t}')
    ##plt.legend()
    #plt.ylabel('<d>')
    #plt.xlabel('step')
    #plt.tight_layout()
    #plt.savefig(f'./difussion_mean.svg', dpi=300, format='svg')
    #plt.savefig(f'./difussion_mean.png', dpi=300, format='png')
    
    #Variance
    #plt.figure(figsize=(8,6))
    #for t in temps:
    for t in [30]:
        L = list(d_results_list[t].keys())
        L.sort()
        y = []
        step_list=[]

        #Compute moments of d
        for step in L:
            if step not in [1,4,10,60]:
                step_list.append(step)
                std = np.array(d_results_list[t][step]).std()
                y.append(std)

        #plt.plot(np.array(step_list), np.array(y), 'o-', c=cmap[t], label=f'{t}')
    ##plt.legend()
    #plt.ylabel(f'{greek_leeterz[18]}')
    #plt.xlabel('step')
    #plt.tight_layout()
    #plt.savefig(f'./difussion_std.svg', dpi=300, format='svg')
    #plt.savefig(f'./difussion_std.png', dpi=300, format='png')

    #Index

    plt.figure(figsize=(8,6))
    #for t in temps:
    for t in [30]:
        L = list(d_results_list[t].keys())
        L.sort()
        y = []
        step_list=[]

        #Compute moments of d
        for step in L:
            if step not in [1,4,10,60]:
                step_list.append(step)
                std = np.array(d_results_list[t][step]).std()
                mean = np.array(d_results_list[t][step]).mean()
                y.append(std**2/mean)
        y = np.array(y)
        plt.plot(np.array(step_list), np.array(y), symbmap[t], c=cmap[t], label=f'{t}')

    plt.plot(np.linspace(0, 61, 100), np.ones(100), '--', c='black', label='D=1')
    plt.yticks(np.linspace(min(y), max(y)*1.1, 5))
    plt.xticks(np.linspace(min(step_list), max(step_list), 5))
    #plt.tick_params(axis='y', labelsize=1.5 * plt.rcParams['font.size'])
    #plt.tick_params(axis='x', labelsize=1.5 * plt.rcParams['font.size'])
    plt.rcParams.update({'font.size': 18})
    plt.gca().yaxis.set_major_formatter(FuncFormatter(format_func))
    ##plt.legend()
    plt.ylabel(f'Index of Disperssion D')
    plt.xlabel('step')
    plt.tight_layout()
    plt.savefig(f'./difussion_index.svg', dpi=300, format='svg')
    plt.savefig(f'./difussion_index.png', dpi=300, format='png')
    
#Activate this to plot the d distributions    
if 0==1:
    d_hist = {t:{} for t in temps}
    
    for t in temps:
        plt.figure(figsize=(8,6))
        
        L = list(d_results_list[t].keys())
        L.sort()

        desired_steps = [2,12,30, 58]
        #desired_steps = [2]
        print(desired_steps)
        
        #Compute moments of d
        for step in L: 
            if step in desired_steps:
                d_hist[t][step] = {}
                d_list = d_results_list[t][step]
                n = len(d_list)


                for d in d_list:
                    try: d_hist[t][step][d] += 1/n
                    except: d_hist[t][step][d] = 1/n

                x = list(d_hist[t][step].keys())
                x.sort()
                x = np.array(x)
                y = np.array([d_hist[t][step][d] for d in x])
                plt.plot(x, y, 'o-', c=cmap[t], marker=symbdict[step], label=f'{step}')


                #Add a Poisson distribution
                x = np.linspace(0,10, 100, dtype=int)
                d_list = np.array(d_list)
                d_list = d_list[d_list<=5]
                #mean = np.array(d_list).mean()
                mean = 0.35
                y = np.array([mean**d * np.e**(-mean)/np.math.factorial(d) for d in x])
                plt.plot(x,y, '--', label=f'Poisson {round(mean,2)}')

        plt.ylim(10**-6, 1.1)

        plt.ylabel('p(d)')
        plt.yscale('log')
        plt.xlabel('d')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f'./pd_{t}.svg', dpi=300, format='svg')
        plt.savefig(f'./pd_{t}.png', dpi=300, format='png')
        plt.xscale('log')
        plt.savefig(f'./pd_{t}_log.svg', dpi=300, format='svg')
        plt.savefig(f'./pd_{t}_log.png', dpi=300, format='png')

     
