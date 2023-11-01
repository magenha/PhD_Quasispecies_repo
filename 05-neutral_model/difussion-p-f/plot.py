"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

L=1000
A=4

#f_values = [1,3,10]
f_values = [3]

f_marker_dict = {1:'o', 3:'^', 10:'+'}

df_sim = pd.read_csv('simulated.csv', header=None)
df_sim.columns=['m']

#Compute the theoretical curve
x = np.arange(0,len(df_sim))
y = (A-1)/A * L * (1-np.exp(-1.0*(A/(A-1))*x/L))

plt.figure();


#Search the string-simulated curve
for f in f_values:
     
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(df['d2'][:5000], f"{f_marker_dict[f]}", markersize=1, label=f'simulated f={f}'); 
    plt.plot(df_sim['m'][:5000], label='Stochastic')
#plt.plot(df_sim['m'],  label='m simulation')
plt.plot(x[:5000],y[:5000],'-', c='black', markersize=5, label='Analytic')


plt.plot(df['d2'][:5000], df['d2'][:5000], label='normal diffusion');
plt.plot(x[:5000],[(A-1)*L/A for i in x][:5000], '--', c='black', label=f'L(A-1)/A'); 

plt.ylabel('<m>')
plt.xlabel('step')
plt.legend(loc='best')
plt.savefig('./example.svg', dpi=300, format='svg')
plt.show()

'''
fig = plt.figure(); 
plt.plot(df_sim['m'], label='m simulation');
plt.plot(x,y, '--', c='black', label='Analytical')
for f in f_values:
    
    ax = fig.gca(); 
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(df['d2'],f"{f_marker_dict[f]}", markersize=1, label=f'simulated f={f}'); 
    #plt.plot(trivialD2, label='alpha=1'); 
    
    ax.set_xscale("log"); 
    ax.set_yscale("log"); 
    plt.ylabel('<d^2>')
    plt.xlabel('t')
    plt.legend(loc='best')



plt.show(); '''
