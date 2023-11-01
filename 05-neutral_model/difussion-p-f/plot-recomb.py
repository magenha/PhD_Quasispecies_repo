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
    x = np.arange(1,len(df)+1)
    plt.plot(x, df['d2'], f"{f_marker_dict[f]}", markersize=3, label=f'recombination f={f}'); 
    
    #plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel('<m>')
    plt.xlabel('t')
    plt.legend(loc='best')
plt.plot(df['d2'], df['d2'], label='alpha=1'); 
plt.plot(x,[(A-1)*L/A for i in x], '--', c='black', label=f'L(A-1)/A')
#plt.yscale('log')
#plt.xscale('log')

plt.savefig('./recomb-fs.png', dpi=300, format='png')
plt.savefig('./recomb-fs.svg', dpi=300, format='svg')
plt.legend(loc='best')
plt.show()