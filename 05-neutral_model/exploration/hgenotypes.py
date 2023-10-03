'''
    This module is intended to contain all tools to generate and mutate individual 
    and population of genotypes. 
''' 
import random
import copy
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import operator

alphabet = ['A', 'C', 'T', 'G']
#alphabet = ['A', 'B', 'C', 'D']

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

protein_wt = ''
our_reference_sequence_qbeta = 'CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG'


our_first_sequence = our_reference_sequence_qbeta
#our_first_sequence = 'AAAAA'
#protein_wt = seq_2_amino(our_first_sequence)[0]

class Genotype():
    
    #Methods
    def __init__(self, fitness, sequence='', id=None):
        '''
            Method to initialize the genotype class
            Input:
                -sequence,  default=''
                -id,     default=None
                -fitness, default=1.0
                
        '''
        
        if len(sequence)==0:
            self.sequence = ''
        else:
            self.sequence = sequence
        if id:
            self.id = id
        else:
            self.id = None
            
        self.L = len(list(self.sequence))
        
        self.fitness = fitness #random.randint(1,100)
        
        
        return None
    def __repr__(self):
        if self.id:
            return f'{self.id} -> {self.sequence}'
        else:
            return f'None -> {self.sequence}'
    def __str__(self):
        if self.id:
            return f'{self.id} -> {self.sequence}'
        else:
            return f'None -> {self.sequence}'
        
    def __copy__(self):
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, number):
            return self.number == other.number
        return False
        
    #
    def mutate(self,n_muts):
        
        global alphabet, protein_wt        
        
        g_new = copy.copy(self)
        
        #Select a (random) locus and mutate using a (random) alphabet letter
        if n_muts:
            locus_list = random.sample(range(len(g_new.sequence)), n_muts)
            for locus in locus_list:
                alph = alphabet.copy()
                alph.remove(g_new.sequence[locus])
                substitution = random.choice(alph)
                g_new.sequence = g_new.sequence[:locus] + substitution + g_new.sequence[locus+1:]
                
        return g_new

#Save an entire population
def save_pop(list_mutants, file_output):
    '''
        File intended to export the population in a collapsed data shape
    '''
    abundance = {}

    for genotype in list_mutants:
        try:
            abundance[genotype.sequence] +=1
        except:
            abundance[genotype.sequence] = 1

    
    with open(file_output, 'w') as f:
        for i, g in enumerate(abundance.keys()):
            abundance_i = abundance[g]
            sequence_i = str(g)
            f.write(f'>{i}-{abundance_i}\n')
            f.write(f'{sequence_i}\n')


