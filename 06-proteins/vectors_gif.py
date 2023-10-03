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
import pickle


#Directory holding all data
dir_path = '/home/samuel/Documents/PhD/Quasispecies/Data/Sequences_filtered/Proteins/'

#Name of the index file
index_name =  dir_path + 'seqs_index.dict'
#Folder contaning the files
mypath = dir_path 

#Output directory
out_dir = './Results/'


#Read index
index_dict = {}

with open(index_name, 'r') as f:
    for line in f:
        L = line.split('\t')
        index_dict[L[1][:-1]] = int(L[0])


with open(out_dir + 'Computations/df.pickle', 'rb') as f:
    df = pickle.load(f)
#print(df)
###########################################################333

temps = [30,43] 
cmap = {30:'b', 43:'r'}

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
t_view,f_view = 0,0
measures = {m:{t:{} for t in temps} for m in ['phi', 'theta']}

#Read th and phi
for t in temps:
    L = list(df[t].columns)
    #print(L)
    L.sort()

    #Need phi & th for certain temp

    with open(f'{out_dir}Computations/phi_{t}', 'rb') as f:
        measures['phi'][t]['normal'] = pickle.load(f)
    with open(f'{out_dir}Computations/phi_{t}_log', 'rb') as f:
        measures['phi'][t]['log'] = pickle.load(f)
    with open(f'{out_dir}Computations/th_{t}', 'rb') as f:
        measures['theta'][t]['normal'] = pickle.load(f)
    with open(f'{out_dir}Computations/th_{t}_log', 'rb') as f:
        measures['theta'][t]['log'] = pickle.load(f)

for aspect in ['normal', 'log']:


    for i in range(1,len(L)):
        step = L[i]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for t in temps:
            phi=measures['phi'][t][aspect]*np.pi/180
            th = measures['theta'][t][aspect]*np.pi/180
            xx = np.sin(th[:i])*np.cos(phi[:i])
            yy = np.sin(th[:i])*np.sin(phi[:i])
            zz = np.cos(th[:i])
            # Create a sphere
            phi_aux , theta_aux = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
            #phi_aux = np.linspace(0, pi/2, 100)
            #theta_aux = np.linspace(0, 2*pi, 100)
            x = r*sin(phi_aux)*cos(theta_aux)
            y = r*sin(phi_aux)*sin(theta_aux)
            z = r*cos(phi_aux)

            #Set colours and render
            ax.scatter(x, y, z, alpha=0.01)

            if t ==temps[0]:
                ax.plot(xx,yy,zz,color="b", alpha=1.0, label=f't={t}')
                ax.scatter(xx[0],yy[0],zz[0],color="b", s=10, alpha=1.0)
                ax.scatter(xx[-1],yy[-1],zz[-1],color="y", s=10, alpha=1.0)
            elif t==temps[1]:
                ax.plot(xx,yy,zz,color="r", alpha=1.0, label=f't={t}')
                ax.scatter(xx[0],yy[0],zz[0],color="b", s=10, alpha=1.0)
                ax.scatter(xx[-1],yy[-1],zz[-1],color="g", s=10, alpha=1.0)

            ax.set_xlim([0,1])
            ax.set_ylim([0,1])
            ax.set_zlim([0,1])
            ax.set_box_aspect((1,1,1))
        
            ax.view_init(t_view, f_view)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'step={step}')
        plt.tight_layout()
        plt.legend()
        plt.savefig(f'{out_dir}Images/Temp/sphere_repr_{t_view}_{f_view}_{step}_merged_{aspect}.png')

    # Create new figure for GIF
    fig, ax = plt.subplots()

    # Adjust figure so GIF does not have extra whitespace
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    ax.axis('off')
    ims = []

    for m in range(1,len(L)):
        step = L[m]
        im = ax.imshow(plt.imread(f'{out_dir}Images/Temp/sphere_repr_{t_view}_{f_view}_{step}_merged_{aspect}.png'), animated = True)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=400)
    ani.save(f'{out_dir}Images/sphere_repr_{t_view}_{f_view}_{aspect}.gif')


#Remove all files in temporal folder