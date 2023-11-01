"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

L=1000
A=4

f_values = ['sqrt', 'linear', 'pow2', 'inverse']
#f_values = ['sqrt']

f_marker_dict = {'sqrt':'o', 3:'^', 'linear':'+', 'pow2':'o', 'inverse':'o'}



plt.figure();

for f in f_values:
     
    df = pd.read_csv(f"dataOut_corrFL_{f}_.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    x = np.arange(1,len(df)+1)
    plt.plot(x,df['d2'], f"{f_marker_dict[f]}", markersize=1, label=f'f {f} m'); 
    
    #plt.plot(df_sim['m'], label='m simulation')
plt.ylabel('<m>')
plt.xlabel('t')
plt.legend(loc='best')
plt.plot(df['d2'], df['d2'], label='alpha=1'); 


plt.savefig('./corrFL.png', dpi=300, format='png')
plt.savefig('./corrFL.svg', dpi=300, format='svg')

plt.show(); 
