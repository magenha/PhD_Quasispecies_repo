"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

L=1000
A=4

f_values = [1,3, 10]

f_marker_dict = {1:'o', 3:'^', 10:'+'}



plt.figure();
#plt.plot(df_sim['m'],  label='m simulation')
#plt.plot(x,y,'--', c='black', label='Analytical')
for f in f_values:
     
    df = pd.read_csv(f"dataOut_f_{f}_recombination.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(df['d2'], f"{f_marker_dict[f]}", markersize=3, label=f'recombination f={f}'); 
    
    #plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel('<d^2>')
    plt.xlabel('t')
    plt.legend(loc='best')
plt.plot(df['d2'], df['d2'], label='alpha=1'); 
#plt.yscale('log')
#plt.xscale('log')


plt.show()