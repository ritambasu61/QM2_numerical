#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 23:09:04 2020

@author: ritambasu
"""
#assuming periodic boundary conditions
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import fixed_quad

hbar = 1; m = 0.01;

# Potential Function
def v(x):
    return (x**4 - 3) * np.exp(- x**2 / 2)



#system parameters
L_start = -5
L_end = 5
L_step = 1
L_arr = [] # to sotre the all L 




# Lists to store the energies for all L an k-values for 1st BZ for all L
Energys_for_all_L = []
K_arr_for_all_L = []
while L_start >= -25 and L_end <= 25:
    # List to store the energies for all K values at a fixed L
    Energys_for_fixed_L = []
    L = (L_end - L_start)
    K_max_1st_BZ = np.pi / L
    K_min_1st_BZ = -np.pi / L
    K_arr = np.linspace(K_min_1st_BZ , K_max_1st_BZ ,25) # for the k values 
    K_arr_for_all_L.append(K_arr)
    for p in range(len(K_arr)):
        # the Dimension of the matrix is dim i.e the no of eigenvalues
        n = 15
        dim = 2*n + 1
        g = 2*np.pi / L
        # Construction of Hamintonian
        H = np.zeros((dim,dim))
        for i in range (dim):
            k = K_arr[p] + (i-n) * g  
            
            for j in range(dim):
                # integrand for the matrix elements
                def integrand(x):
                    return v(x) * np.cos((i-j)*g*x)

                if i == j:
                    H[i][j] = (hbar * k)**2 / (2*m)
                elif j > i:
                       H[i][j] = (1/L) * fixed_quad(integrand , -L/2 , L/2)[0]
        for i in range(dim):
            for j in range(dim):
                if i > j:
                    H[i][j] = H[j][i]     
        E = np.linalg.eigh(H)[0]
        Energys_for_fixed_L.append(E)
    
    Energys_for_all_L.append(Energys_for_fixed_L)
    L_arr.append(L)
    L_start = L_start - L_step/2
    L_end = L_end + L_step/2





# Band Structure for L=10
En_L_10 = np.transpose(np.array(Energys_for_all_L[0]))
for i in range(6):     
    plt.plot(K_arr_for_all_L[0] , En_L_10[i])
print('Band Gap Near the 1st BZ for m=',m,' : ' ,min(En_L_10[1]) - max(En_L_10[0]) )
plt.xlabel('k' , fontsize = 16 , color = 'Black')
plt.ylabel('E(k)', fontsize = 16, color = 'Black')
plt.title('Band Structure for L=10', fontsize = 18 , color = 'Black')
plt.show()




#plotting the minimum of each band as function of L
# store the minimum of each band for all L
min_of_bands = []
l1 = len(Energys_for_all_L)
l2 = len(E)
for i in range(l1):
    # store the minimum of each band for a fixed L
    min_bands_fixed_L = []
    for j in range(l2):
        min_bands_fixed_L.append(min(np.transpose(np.array(Energys_for_all_L[i]))[j]))
    min_of_bands.append(min_bands_fixed_L)

#plotting minimum value of first 30 bands as a function of L
num = 30
# Plotting the minimum of even bands 
for i in range(int(num/2)):     
    plt.plot(L_arr , np.transpose(np.array(min_of_bands))[2*i],label='n = {0:d}' .format((2*i)+1))
plt.xlabel('L' , fontsize = 16, color = 'Black')
plt.ylabel('E(L)', fontsize = 16 , color = 'Black')
plt.title('Minimum of Band Energy Values E(L) vs L\n for even bands', fontsize = 18 , color = 'Black')
plt.grid(color='Black', linestyle=':', linewidth=0.5)
plt.legend(fontsize=7,loc='best')
plt.show()

# Plotting the minimum of odd bands 
for i in range(int(num/2)):     
    plt.plot(L_arr , np.transpose(np.array(min_of_bands))[2*i+1],label='n = {0:d}' .format((2*i)+1))
plt.xlabel('L' , fontsize = 16 , color = 'Black')
plt.ylabel('E(L)', fontsize = 16 , color = 'Black')
plt.title('Minimum of Band Energy Values E(L) vs L\n for odd bands', fontsize = 18 , color = 'Black')
plt.grid(color='Black', linestyle=':', linewidth=0.5)
plt.legend(fontsize=7,loc='best')
plt.show()
   