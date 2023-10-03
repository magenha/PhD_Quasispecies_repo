'''
    This code simulates the trajectory in Genotype Space by mutating a sequence.
'''

import numpy as np
import random
import hgenotypes as hg
import matplotlib.pyplot as plt

def distance(seq, seq_ref):
    '''
        Function that returns the Humming distance between
        two sequences
        Input:
            -sequence to compare
            -reference sequence
        Return:
            -int : distance
    '''
    d=0
    for char in range(len(seq_ref)):
        if seq[char] != seq_ref[char]:
            d+=1
    return d

def do_mutation(sequence):
    '''
        Function that, given a sequence, returns a mutation
    '''
    global alphabet
    L = len(sequence)
    i = random.randint(0,L-1)
    nucl = random.choice(alphabet)
    while nucl == sequence[i]:
        nucl = random.choice(alphabet)
    if i<L:
        seq_new = sequence[:i] + nucl + sequence[i+1:]
    else:
        seq_new = sequence[:L-1] + nucl
    return seq_new

def calculateParametersLSR(x, y):
     
	"""
	Function that calculates w* and b* for a least square regression and the standard error for b*
	Inputs:
		x: array containing the x-values in log-scale
		y: array containing the y-values in log-scale
	Outputs:
		b: float containing the intercept term of the regression line
		w: float containing the slope of the regression line
		SEw: float containing the standard error of b*
	"""

	meanTarget = np.mean(y)
	meanFeature = np.mean(x)

	centeredTarget = y - meanTarget
	centeredFeature = x - meanFeature

	w = (centeredFeature @ centeredTarget)/(centeredFeature @ centeredFeature)

	b = meanTarget - w *  meanFeature

	# Standard error
	yHat = b + w * x
	n = x.shape[0]

	SEw = np.sqrt((1/(n - 2)) * ((np.sum(np.power((y - yHat), 2)))/(np.sum(np.power((x - meanFeature), 2)))))



	return b, w, SEw



#Initial sequence
#First term in hg.Genotype() is the fitness (number of offspring when replicated)
#Second term in hg.Genotype() function is the given sequence. 
L = 10**2
initial_sequence = ''.join(['A']*L)
s0 = hg.Genotype(1, initial_sequence)
alphabet = ['A', 'C', 'G', 'T']
t_max = 10**3
n_iterations = 10**3
mu = 1*10**-3      #mutation rate -> prob. of mutation per locus
greek_leeterz = [chr(code) for code in range(945,970)]


f_dict = {s0.sequence: 1}  #Initial fitness dictionary  {sequence: f}
p = 0.8   #Fraction of null-fitness points in Genotype Space. (measure of the holey in space)

#This is for further use
time = np.array(range(t_max)) 
distance_t = {t:[] for t in range(t_max)}


print('Initiating rounds')
for i in range(n_iterations):  
    print(i)
    #Define the initial sequence -> equals to hg.qBeta  
    s = s0
    distance_t[0].append(distance(s.sequence, s0.sequence))


    for t in range(1,t_max):
        #Define the number of mutations that the sequence will experiment.
        #This can be set = 1
        jump = np.random.poisson(mu*L)


        #Check if the sequence is able to replicate or mutate.
        if f_dict[s.sequence]:
            #Mutate the sequence
            #method mutate in hg library selects n_muts loci of the sequence and changes the letter
            s = s.mutate(n_muts=jump)

            #This is  variant to use the function do_mutate, that realizes a single mutation
            #s = hg.Genotype(1, do_mutation(s.sequence))  

            #Check if the new sequence is in fitness dict
            if s.sequence in f_dict.keys():
                pass
            else:
                #Add the new sequence to dict, and assign value of !0 fitness with probability p
                r = random.random()
                if r<=p:
                    f_dict[s.sequence] = 1
                else:
                    f_dict[s.sequence] = 0
                
        distance_t[t].append(distance(s.sequence, s0.sequence))


position = [np.array(distance_t[t]).mean() for t in range(t_max)]
variance = [(np.array(distance_t[t]).std())**2 for t in range(t_max)]
plt.plot(time, position, 'o-', label=f'{greek_leeterz[11]}={mu}')
#plt.plot(time, variance)

t_fit_max = t_max

b,w,E = calculateParametersLSR(np.log2(np.array(time[1:t_fit_max])), np.log2(np.array(position[1:t_fit_max])))
x_aster = time[1:t_fit_max]
y_aster = 2**b*x_aster**(w)
w = round(w,3)
E = round(E,3)
theor_alpha = 1.90
#plt.plot(x_aster, y_aster, '--', label=f'alpha={w}+-{E}')
#plt.plot(x_aster, 2**b*x_aster**(1.94), '--', label=f'alpha={1.9} (for p=1)')
plt.plot(x_aster, 2**b*x_aster**(1.0), '--', label=f'alpha={1.0} (normal diffusion)')
plt.ylabel('<d^2>')
plt.xlabel('t')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('./d_2_t_fitted.png', dpi=300, format='png')