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
index_name = 'seqs_index.dict'

index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])
        
print('Loading data', datetime.now())
WD = os.getcwd()
mypath = WD + '/data'
file_output = WD + '/'+ 'seqs_index.dict'

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
    
    file_name = WD + '/data/' + a
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


temps = [30, 43]

df = {'': {'': {}, 'norm-euc': {}, 'norm-1': {}},
      'log': {'': {}, 'norm-euc': {}, 'norm-1': {}}}
        
print('Adding normalizations')
#Upload data to DataFrame
for t in temps:

    df[''][''][t] = pd.DataFrame([])
    df[''][''][t] = pd.DataFrame.from_dict(data_dict[t])
    df[''][''][t] = df[''][''][t].fillna(0)
    
    #Create log unnormalized
    df_t = df[''][''][t]
    df_2 = df_t + 1
    df['log'][''][t] = df_2.apply('log')
    
    
    #Sort columns of abundance and log
    L = list(df_t.columns)
    L.sort()
    
    #Create empty containers
    for k_1 in df.keys():
        for k_2 in ['norm-euc', 'norm-1']:
            df[k_1][k_2][t] = pd.DataFrame([])
            df[k_1][k_2][t].index = df_t.index
    
    for x in L:
        for k_1 in df.keys():
            arr = np.array(df[k_1][''][t][x].tolist())
            normalized_arr = arr / np.sqrt(np.sum(arr**2))
            relative_arr = arr / np.sum(arr)
            df[k_1]['norm-euc'][t] = normalized_arr
            df[k_1]['norm-1'][t] = normalized_arr

def jaccknife_df(df_in, boxes_number = 10):
    '''
    Function that splits the df in a number of given boxes
    (default=10). 
    
    For the process it desamples the collapsed column, and 
    randomly builds the boxes.
    
    Output: list of the dfs built using boxes_number-1
            (len(list)=boxes_number)
    '''
    column = df_in.to_dict()
    index = df_in.index
    
    desampled = []

    for element in column:
        desampled.append([element for i in range(int(column[element]))])
    desampled = list(itertools.chain(*desampled))
    
    
    #Randomly shuffle the desampled list
    desampled_random = []
    random.seed(datetime.now().timestamp())
    for i in range(len(desampled)):
        j = random.randint(0, len(desampled)-1)
        desampled_random.append(desampled[j])
        desampled.pop(j)   
    desampled = desampled_random

    boxes = {}
    boxes_df = {}
    
    N = len(desampled)
    box_size = int(N//boxes_number)

    
    # Create boxes
    for i in range(boxes_number):
        boxes[i] = desampled[i*box_size:(i+1)*box_size]

 
    for i in range(boxes_number):
        df_x = pd.DataFrame([], columns = ['seq'])
        for j in range(boxes_number):
            if j != i:
                df_x = pd.concat([df_x, pd.DataFrame(boxes[j], columns= ['seq'])])
           
        df_x = df_x.pivot_table(columns=['seq'], aggfunc='size')
        df_aux = pd.DataFrame([], index=index, columns=['seq']) 
        df_aux = pd.concat([df_x,df_aux], axis=0)
        df_aux = df_aux[~df_aux.index.duplicated(keep='first')]
        df_aux = df_aux.sort_index()
        df_aux = df_aux[0].fillna(0)
        boxes_df[i] = df_aux

        #print(df_aux[0])
        
    return list(boxes_df.values())


print('Starting computation of phi and theta with Jackknifing')
#Compute phi: consecutive steps angle
realization = ''
#for t in temps:
#Select temperature
for t in temps:

    L = list(df[''][''][t].columns)
    L.sort()


    populations = []

    #for step in steps
    for j, step in enumerate(L):
        print(f'Jackknifing {realization} temp {t} step {step}')
        populations.append(jaccknife_df(df[''][''][t][step]))

    results_phi = np.zeros((len(L)-1,2))
    results_th = np.zeros((len(L),2))


    for j,step in enumerate(L):
        if (j < len(L)-1):
            cos_phi = []
            cos_th = []
            for i in range(len(populations[j])):
                p = populations[j][i]
                p_n = populations[j+1][i]
                
                a = np.array(p.tolist())
                b = np.array(p_n.tolist())
                c = np.array(populations[0][i])

                cos_phi_aux = np.dot(a,b)/(np.sqrt(np.sum(a**2))*np.sqrt(np.sum(b**2)))
                cos_phi.append(cos_phi_aux)

                cos_th_aux = np.dot(a,c)/(np.sqrt(np.sum(a**2))*np.sqrt(np.sum(c**2)))
                cos_th.append(cos_th_aux)
            N = len(cos_phi)
            cos_phi = np.array(cos_phi)
            step_norm = L[j+1]-L[j]
            cos_phi = np.around(cos_phi, decimals=6)
            phi = np.arccos(np.array(cos_phi))*180/np.pi/step_norm
            phi_mean_quantity = np.mean(phi)

            cos_th = np.array(cos_th)
            cos_th = np.around(cos_th, decimals=6)
            th = np.arccos(np.array(cos_th))*180/np.pi
            th_mean_quantity = np.mean(th)

            results_phi[j][0] = phi_mean_quantity
            results_phi[j][1] = len(populations) * np.sqrt(np.sum(np.power((phi-phi_mean_quantity),2))/(N-1))
            results_th[j][0] = th_mean_quantity
            results_th[j][1] = len(populations) * np.sqrt(np.sum(np.power((th-th_mean_quantity),2))/(N-1))
            
            #Try to see error
            if results_phi[j][1] >=5.0 :
                with open(f'./error_{t}_{j}.txt', 'w') as f:
                    for k, element in enumerate(cos_phi):
                        f.write(f'{cos_phi[k]}, {phi[k]}\n')
                    f.write(f'mean={phi_mean_quantity}, std={results_phi[j][1]}')
                    
                    f.write(f'{phi-phi_mean_quantity}')
                    f.write(f'{np.power(phi-phi_mean_quantity,2)}')
                    f.write(f'{np.power(phi-phi_mean_quantity,2)/(N-1)}')
                    f.write(f'{np.sum(np.power((phi-phi_mean_quantity),2))/(N-1)}')
                    f.write(f'{len(populations) *np.sqrt(np.sum(np.power((phi-phi_mean_quantity),2))/(N-1))}')
        else:
            cos_th = []
            for i in range(len(populations[j])):
                p = populations[j][i]
                a = np.array(p.tolist())
                c = np.array(populations[0][i])

                cos_th_aux = np.dot(a,c)/(np.sqrt(np.sum(a**2))*np.sqrt(np.sum(c**2)))
                cos_th.append(cos_th_aux)

            cos_th = np.array(cos_th)
            cos_th = np.around(cos_th, decimals=4)
            th = np.arccos(np.array(cos_th))*180/np.pi
            th_mean_quantity = np.mean(th)
            th_std_quantity = len(populations) * np.std(th)
            results_th[j][0] = th_mean_quantity
            results_th[j][1] = th_std_quantity    
    #Export the results array to files
    with open(f'./computes_results/th_{t}_{realization}.np', 'wb') as f:
        np.save(f, results_th)

    with open(f'./computes_results/phi_{t}_{realization}.np', 'wb') as f:
        np.save(f, results_phi)

    #Save steps
    with open(f'./computes_results/steps_phi_{t}_{realization}.np', 'wb') as f:
        np.save(f,np.array(L[:-1]))
    with open(f'./computes_results/steps_theta_{t}_{realization}.np', 'wb') as f:
        np.save(f,np.array(L))


print('Building figures and plots')

#Load data
quantities = {'phi' : {},
              'th' : {}}
temps = [30,43]

cmap = {30:'b', 43:'r'}
steps = {q:{t:[] for t in temps} for q in quantities.keys()}

for q in list(quantities.keys()):
    for t in temps:  
        with open(f'./computes_results/{q}_{t}_{realization}.np', 'rb') as f:
            quantities[q][t] = np.load(f)
        with open(f'./computes_results/steps_{q}_{t}_{realization}.np', 'rb') as f:
            steps[q][t] = np.load(f)
                

greek_letterz=[chr(code) for code in range(945,970)]           
letter_theta = greek_letterz[7]
letter_phi = greek_letterz[21]

#Initiate fig
#Create plot for phi, and all temps
#Save plot

plt.figure()

for t in temps:
    x = steps['phi'][t][:,0]
    y = quantities['phi'][t][:,0]
    y_err = quantities['phi'][t][:,1]
    

    plt.plot(x,y, 'o-', label=f't={t}', c=cmap[t])
    plt.errorbar(x,y,y_err)
    plt.ylabel(f'{letter_phi}')
    plt.xlabel('step')
    plt.yscale('linear')
    plt.legend()
    plt.tight_layout()
    plt.title(f'phi')
plt.savefig(f'./pics_results/phi_{realization}.jpg')

plt.figure()

for t in temps:
    x = steps['phi'][t][:,0]
    y_aux = quantities['phi'][t][:,0]
    y = [y_aux[0]]
    for i in range(1,len(y_aux)):
        y.append(y[i-1]+y_aux[i])


    #y_err = quantities['phi'][t][:,1]
    

    plt.plot(x,y, 'o-', label=f't={t}', c=cmap[t])
    #plt.errorbar(x,y,y_err)
    plt.ylabel(f'sum({letter_phi})')
    plt.xlabel('step')
    plt.yscale('linear')
    plt.legend()
    plt.tight_layout()
    plt.title(f'phi')
plt.savefig(f'./pics_results/phi_cummulative.jpg')


#Initiate fig
#Create plot for theta, and all temps
#Save plot

plt.figure(figsize=(12,9))

for t in temps:
    x = quantities['th'][t][:,0]
    y = quantities['th'][t][:,0]
    y_err = quantities['th'][t][:,1]
    

    plt.plot(x,y, '-r')
    plt.scatter(x,y, label=f't={t}')
    plt.errorbar(x,y,y_err)
    plt.ylabel(f'{letter_theta}')
    plt.xlabel('step')
    plt.yscale('linear')
    plt.legend()  
    plt.tight_layout()
    plt.title(f'theta')
plt.savefig(f'./pics_results/th_{realization}.jpg')
