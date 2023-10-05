"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

L=1000
A=4

f_values = [1,3,10]

f_marker_dict = {1:'o', 3:'^', 10:'+'}

df_sim = pd.read_csv('simulated.csv', header=None)
df_sim.columns=['m']

plt.figure()
for f in f_values:
    
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']

    dR = df['d2'].tolist()
    tMax = len(dR); 
    trivialD2 = [tt for tt in range(tMax)]; 

    plt.plot(df['variance'].tolist(), f"{f_marker_dict[f]}", markersize=3, label=f'f={f}')
    plt.ylabel('Variance')
    plt.xlabel('t')
    plt.legend(loc='best')

plt.figure()
for f in f_values:
    
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']

    plt.plot(df['mean'].tolist(), f"{f_marker_dict[f]}", markersize=3, label=f'f={f}')
    plt.ylabel('<d>')
    plt.xlabel('t')
    plt.legend(loc='best')



x = np.arange(0,len(df_sim))
y = (A-1)/A * L * (1-np.exp(-1.0*(A/(A-1))*x/L))

plt.figure();
plt.plot(df_sim['m'],  label='m simulation')
plt.plot(x,y,'--', c='black', label='Analytical')
for f in f_values:
     
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(dR, f"{f_marker_dict[f]}", markersize=3, label=f'simulated f={f}'); 
    #plt.plot(trivialD2, label='alpha=1'); 
    #plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel('<d^2>')
    plt.xlabel('t')
    plt.legend(loc='best')

fig = plt.figure(); 
plt.plot(df_sim['m'], label='m simulation');
plt.plot(x,y, '--', c='black', label='Analytical')
for f in f_values:
    
    ax = fig.gca(); 
    df = pd.read_csv(f"dataOut_f_{f}.csv", header=None)
    df.columns = ['mean', 'variance', 'd2']
    plt.plot(dR,f"{f_marker_dict[f]}", markersize=3, label=f'simulated f={f}'); 
    #plt.plot(trivialD2, label='alpha=1'); 
    
    ax.set_xscale("log"); 
    ax.set_yscale("log"); 
    plt.ylabel('<d^2>')
    plt.xlabel('t')
    plt.legend(loc='best')



plt.show(); 
