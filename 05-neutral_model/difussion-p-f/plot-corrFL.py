"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

L=1000
A=4

f_values = ['sqrt', 'linear', 'pow2_', 'inverse_']
#f_values = ['sqrt']

f_marker_dict = {'sqrt':'o', 3:'^', 'linear':'+', 'pow2_':'o', 'inverse_':'o'}



plt.figure();

for f in f_values:
     
    df = pd.read_csv(f"dataOut_corrFL_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(df['d2'], f"{f_marker_dict[f]}", markersize=1, label=f'f {f} m'); 
    
    #plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel('<d^2>')
    plt.xlabel('t')
    plt.legend(loc='best')
plt.plot(df['d2'], df['d2'], label='alpha=1'); 

plt.yscale('log')
plt.xscale('log')
plt.show(); 
