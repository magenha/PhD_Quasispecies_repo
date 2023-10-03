"""plot.py"""


import numpy as np; 
import pandas as pd;
import matplotlib.pyplot as plt; 
import random as rnd; 
from copy import copy; 

f_value = 10

df = pd.read_csv(f"dataOut_f_{f_value}.csv", header=None)
df.columns = ['mean', 'variance', 'd2']

dR = df['d2'].tolist()
tMax = len(dR); 
trivialD2 = [tt for tt in range(tMax)]; 


plt.figure()
plt.plot(df['variance'].tolist())
plt.ylabel('Variance')
plt.xlabel('t')

plt.figure()
plt.plot(df['mean'].tolist())
plt.ylabel('<d>')
plt.xlabel('t')

plt.figure(); 
plt.plot(dR, label='simulated'); 
plt.plot(trivialD2, label='alpha=1'); 
plt.ylabel('<d^2>')
plt.xlabel('t')
plt.legend(loc='best')

fig = plt.figure(); 
ax = fig.gca(); 
plt.plot(dR,label='simulated'); 
plt.plot(trivialD2, label='alpha=1'); 
ax.set_xscale("log"); 
ax.set_yscale("log"); 
plt.ylabel('<d^2>')
plt.xlabel('t')
plt.legend(loc='best')



plt.show(); 
