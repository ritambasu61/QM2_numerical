#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 14:32:28 2020

@author: ritambasu
"""

import numpy as np
import matplotlib.pyplot as plt

#defining the potential energy
def V(x):
    return((x**4-3)*np.exp(-0.5*x**2))


def hamiltonian(V, t0, tn, n=400):
    x= np.linspace(t0, tn, n+2)
    dx = x[1] - x[0]
    kinetic= np.identity(n)
    T = np.identity(n+1)#kinetic energy matrix
    kinetic= (-2*kinetic) + T[1:n+1, 0:n] + T[0:n, 1:n+1]
    kinetic= (-0.5 *kinetic)/dx**2
    V_diag= np.diag(V(x[1:n+1]))#potential after diagonalizing it
    M =  kinetic + V_diag#total diagonalize hamiltonian(kinetic+potential)
    eigenvalue = np.linalg.eigvalsh(M)#solving it,s eigenvalues
    eigenvalue = np.sort(eigenvalue)#sorting eigenvalus in increasing order
    return(eigenvalue)


# Defining Problem paramaters
e1 = hamiltonian(V, -5, 5)#eigenvalues for given potential
print("First 5 Energy levels = ", e1[0:5])#first 5 energy levels within a box of size L = 10 in between -5 to 5
t = np.linspace(5, 25, int((25-5)/0.1)+1)
e = []
for i in t:
    ei = hamiltonian(V, -1*i, i)
    e.append(ei[0:20])

#plotting of even states
e = np.array(e)
for i in range(10):
    plt.plot(2*t, e[:, 2*i], label='n = {0:d}' .format(2*i))
plt.ylim([-1, 5])
plt.title("Even States")
plt.xlabel("$k$",fontsize=14)
plt.ylabel("$E_n$",fontsize=14)
plt.legend(fontsize=8,loc='best')
plt.show()

#plotting of odd states
for i in range(9):
    plt.plot(2*t, e[:, (2*i)+1], label='n = {0:d}' .format((2*i)+1))
plt.ylim([-2.5, 5])
plt.title("Odd States")
plt.xlabel("$k$",fontsize=14)
plt.ylabel("$E_n$",fontsize=14)
plt.legend(fontsize=7,loc='best')
plt.show()
