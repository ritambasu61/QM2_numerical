#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun7  22:16:46 2020

@author: ritambasu
"""

import numpy as np 
import matplotlib.pyplot as plt 

#===============================================================================
# Define parameters:
I = 300 
hcut = 1 
x_min = 0
x_max =  2*np.pi
N = 801
# potential:
def v(x):
	return 1-np.cos(x)

#===============================================================================
x = np.linspace(x_min,x_max,N)
h = x[1]-x[0]
factor = (hcut**2)/(2*I*(h**2)) 

#defining kinetic energy and potential enrgy matrix
T = np.zeros((N-1, N-1))
V = np.zeros((N-1, N-1))
T[0][N-2] = 1
T[N-2][0] = 1 

for i in range(len(T)):
	for j in range(len(T)):
		if (i == j):
			T[i][j] = -2
			V[i][j] = v(x[i+1])
		elif ( abs(i-j) == 1):
			T[i][j] = 1
#defining the hamiltonian
H = - factor*T + V 

#solving it's eigenvalues
EV,Evec = np.linalg.eig(H)
EV  = np.sort(EV) #sorting the eigenvalues          
state_no = np.arange(0,300,1)
print('The energies of ground st:',EV[0],'\nThe energy of 50th state:', EV[49],'\nThe energy of 200th state:',EV[199])

#plotting the energy eigenvalues
plt.plot(state_no,EV[:300],'r')
plt.grid(color='b', linestyle=':', linewidth=1.0)
plt.xlabel('State number(n)',size=14)
plt.ylabel('Energy(E(n))',size=14)
plt.title('Energy vs state number', size=18)
plt.show()



#Finding 'From which state the eigenvalues are very closed' i.e degenerate
count=0
e_degenerate=0
for i in range(len(EV)):
    if abs(EV[i+1]-EV[i])<=0.01:#with respect to machine precision
        count=1
        break
if count==1:
    e_degenerate=EV[i]
    print('Degenaracy start from ',i,'th state with energy:',e_degenerate)
else:
    print('No degeneracy found')
        

    

#plottings
x=np.arange(0,101,0.1)
degen_energy=(np.zeros(len(x))+1)*e_degenerate        
plt.scatter(state_no[:100],EV[:100],s=2.0)
plt.plot(x,degen_energy,color='Red',label='degeneracy start')
plt.grid(color='Black',axis='x',linestyle=':', linewidth=1.0)
plt.xlabel('State number(n)',size=14)
plt.ylabel('Energy(E(n))',size=14)
plt.title('Energy vs state number', size=18)
plt.legend(loc='best')
plt.show()
