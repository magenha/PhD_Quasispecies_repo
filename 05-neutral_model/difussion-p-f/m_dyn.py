
import random
import numpy as np
import matplotlib.pyplot as plt

alphabet = ['A', 'C', 'G', 'T']

L=1000

t_max = 10000
A = len(alphabet)
d_H = np.zeros(t_max)


for i in range(t_max-1):
    m = d_H[i]
    p_p = 1-m/L
    p_m = (1/(A-1))*m/L 
    p_s = 1- p_p - p_m

    r = random.random()

    if r<p_p:
        m+=1
    elif r<p_p+p_m and m>0:
        m+=-1
    else:
        pass
    d_H[i+1] = m


with open(f'./simulated.csv', 'w') as f:
    for i in range(len(d_H)):
        f.write(f'{d_H[i]}\n')


plt.figure()
plt.plot(range(t_max), d_H)
plt.plot(range(t_max), [L*(A-1)/A for i in range(t_max)], '--', label='K=L(A-1)/A')
plt.plot(range(t_max), [0.5*L*(A-1)/A for i in range(t_max)], '--', label='K/2')
plt.plot(range(t_max), [L*(A-1)/(2*A-3) for i in range(t_max)], '--', label='p+=p0')
plt.xlabel('t')
plt.ylabel('m(t)')
#plt.yscale('log')
#plt.xscale('log')

plt.legend(loc='best')

plt.figure()

d_t = 200

d_d_H = [(d_H[i:i+d_t//2].mean()-d_H[i-d_t//2:i].mean())/(d_t) for i in range(d_t//2,t_max-d_t//2)]
 
plt.plot(range(d_t//2,t_max-d_t//2), d_d_H, label=f'D(d_H)/dt, with dt={d_t}')
#plt.plot(range(d_t//2,t_max-d_t//2), [1-A/(A-1)*r/L for r in range(d_t//2,t_max-d_t//2)], label=f'theor')
plt.legend(loc='best')
plt.ylabel('d/dt d_H')
plt.xlabel('t')
plt.xscale('log')
'''
plt.figure()
for i in range(t_max-1):
    if d_H[i]>d_H[i+1]:
        plt.plot([d_H[i]/L, d_H[i]/L], [d_H[i]/L, d_H[i+1]/L], '-', c='red')
    elif d_H[i]<d_H[i+1]:
        plt.plot([d_H[i]/L, d_H[i]/L], [d_H[i]/L, d_H[i+1]/L], '-', c='blue')
    else:
        plt.plot([d_H[i]/L, d_H[i]/L], [d_H[i]/L, d_H[i+1]/L], '-', c='black')
    plt.plot(d_H[i]/L, d_H[i+1]/L, 'o', c='black')
    plt.plot([d_H[i]/L, d_H[i+1]/L], [d_H[i+1]/L, d_H[i+1]/L], '-', c='green')
    #plt.plot([d_H[i]/L, d_H[i+1]/L], [d_H[i+1]/L, d_H[i+1]/L], '-', c='black')
plt.plot(d_H/L, d_H/L, '--')
plt.ylabel('m(t+1)')
plt.xlabel('m(t)')

'''
plt.show()