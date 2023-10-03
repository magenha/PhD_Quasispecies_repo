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
import os
from os.path import isfile, join
from operator import itemgetter
import operator


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


seqs_list = {}
for ikey in index_dict.keys():
    seqs_list[index_dict[ikey]] = ikey

codon_2_seq = {
    'UUU': 'F', 
    'UUC': 'F',
    'UUA': 'L',
    'UUG': 'L',
    'UCU': 'S', 
    'UCC': 'S',
    'UCA': 'S',
    'UCG': 'S',
    'UAU': 'Y', 
    'UAC': 'Y',
    'UAA': 'Stop',
    'UAG': 'Stop',
    'UGU': 'C', 
    'UGC': 'C',
    'UGA': 'Stop',
    'UGG': 'W',
    'CUU': 'L', 
    'CUC': 'L',
    'CUA': 'L',
    'CUG': 'L',
    'CCU': 'P', 
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'CAU': 'H', 
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGU': 'R', 
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AUU': 'I', 
    'AUC': 'I',
    'AUA': 'I',
    'AUG': 'M',
    'ACU': 'T', 
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'AAU': 'N', 
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGU': 'S', 
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GUU': 'V', 
    'GUC': 'V',
    'GUA': 'V',
    'GUG': 'V',
    'GCU': 'A', 
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'GAU': 'D', 
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'GGU': 'G', 
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G'
}
def seq_2_amino(seq, initial=0):
    '''
    Function that given a sequence, return the aminoacid that it expresses.
    
    Input:
        -sequence of DNA or RNA 
        -initial position of reading (ORF) (Optional)
    Output:
        -list of aminoacids
    '''
    x = seq#[::-1]
    x = str(x[initial:]).upper()
    x = x.replace('T', 'U')
    x = list(x)
    aminos = []
    for i in range(len(x)//3):
        codon=''.join(x[3*i:3*i+3])
        codon = str(codon)
        try:
            amin = codon_2_seq[codon]
            #print(f'aminoacid={amin}')
            aminos.append(amin)
        except:
            print(f'codon {codon} not in dict')
    aminos = ''.join(aminos)
    aminos = aminos.split('Stop')
    return aminos

def compare_proteins(amin_1, amin_2):
    '''
    Function that compares two aminoacid and returns the differences.
    Input:
        -amin_1 & amin_2, aminoacids to compare
    Output:
        -dict having position: nucl_amin_1, nucl_amin_2
    '''
    results = {}
    
    if len(amin_1) != len(amin_2):
        #print('Aminoacids of different length!')
        return None
    else:
        for i in range(len(amin_1)):
            if amin_1[i] != amin_2[i]:
                results[i] = [amin_2[i], amin_1[i]]
    return results

def t_step_toarx(t,step):
    '''
        Function that given a temperature and step,
        returns the correct archive name
    '''
    
    if t==30:
        if step<10:
            return f'30-1-0{step}-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed'
        else:
            return f'30-1-{step}-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed'
    elif t==43:
        if step<10:
            return f'c43-p0{step}-3-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed'
        else:
            return f'c43-p{step}-3-EL3_all_all_trim_merged_filter_sort_filter_length_collapsed'
    else:
        return None
    
temps = [30,43]
#amino_ab_evolution = {t:{am: [] for am in list(my_top_aminos.keys())[:11]} for t in temps}
for t in temps:
    
    df_t = df[t]
    L = list(df_t.columns)
    L.sort()
    
    for step in L:
        print(t,step)
        amino_dict = {}         #contains amino: abundance
        amino_mapping = {}      #contains amino: [seqs]
        seqs = list(df_t[step].index)
        seqs = [seqs_list[s] for s in seqs]
        
        print('Looking for aminoacids')
        for sequence in seqs[:100]:
            amino = seq_2_amino(sequence)[0]

            if amino not in list(amino_dict.keys()):
                amino_dict[amino] = df_t[step].loc[index_dict[sequence]]
                amino_mapping[amino] = [sequence]
            else:
                amino_dict[amino] +=  df_t[step].loc[index_dict[sequence]]
                if sequence not in amino_mapping[amino]:
                    amino_mapping[amino].append(sequence)
        print('Saving mapping and aminoacid abundance')
        output_file = '/home/samuel/Documents/PhD/Quasispecies/Data/Sequences_filtered/aminoacids/'
        output_file1 = output_file+t_step_toarx(t,step)
        output_file2 = output_file+'mapping_'+ t_step_toarx(t,step)
        
        with open(output_file1, 'w') as f:    #write the abundances of aminoacids
            for i,k in enumerate(amino_dict.keys()):
                ab = int(amino_dict[k])
                f.write(f'>{i}-{ab}\n')
                f.write(k+'\n')
        
        with open(output_file2, 'w') as f:    #write the abundances of aminoacids
            for i,k in enumerate(amino_mapping.keys()):
                rna_list = [index_dict[s] for s in amino_mapping[k]]
                f.write(f'{k}:{rna_list}\n')
        
