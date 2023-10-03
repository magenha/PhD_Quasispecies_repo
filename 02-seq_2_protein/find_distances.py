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


seqs_list = list(index_dict.keys())

b = seqs_list[0]
c = seqs_list[1]
a = seqs_list[2]
x = seqs_list[12]

d_ab=0
d_ac=0
d_bc=0
d_ax = 0
d_bx = 0
d_cx = 0
for char in range(len(a)):
    if a[char] != b[char]:
        print(f'Diference a -> b: {char+1} {a[char]} {b[char]}')
        d_ab+=1
    if a[char] != c[char]:
        print(f'Diference a -> c: {char+1} {a[char]} {c[char]}')
        d_ac+=1
    if b[char] != c[char]:
        print(f'Diference c -> b: {char+1} {c[char]} {b[char]}')
        d_bc+=1
    if b[char] != x[char]:
        print(f'Diference b -> x: {char+1} {b[char]} {x[char]}')
        d_bx+=1
    if a[char] != x[char]:
        print(f'Diference a -> x: {char+1} {a[char]} {x[char]}')
        d_ax+=1

print(f'a->b = {d_ab}')
print(f'a->c = {d_ac}')
print(f'b->c = {d_bc}')
print(f'b->x = {d_bx}')
print(f'a->x = {d_ax}')

seq_x = 'CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCATAAAGCGTTGAAACTTTG'

seq_x = a[:28] + 'G' + a[29:]

for n,s in enumerate(seqs_list):
    if s==seq_x:
        print('The missing sequence between A and C is', n, 'in our dataset')

