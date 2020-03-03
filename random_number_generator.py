#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:31:09 2019

@author: florian
"""

import numpy as np
import matplotlib.pyplot as plt


try:
    plt.close('all')
except:
    pass

outdir = "Plots/"


a, b, c= 11, 0, 2**2

#a, b, c= 65539, 0, 2**31


n=100000
X   = np.zeros(n)
U   = np.zeros(n)
N   = np.zeros(n)
X[0] = 0.5
U[0] = X[0]/c


for i in range(n-1):
    X[i+1] = np.mod((a*X[i]+b),c)
    U[i+1] = X[i+1]/c
    
for i in range(2,n):
    N[i] = (-2*np.log(U[i-1]))**(1/2)*np.cos(2*np.pi*U[i-2])
    
print(U)
plt.figure(1)
plt.hist(U)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Value',fontsize=14)
plt.ylabel('Frequency',fontsize=14)
plt.title('Histogram of U',fontsize=16)
plt.savefig(outdir + "Random_number_generator_hist_U_a_"+str(a)+".png", dpi = 300)

plt.figure(2)
plt.hist(N,bins=100)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Value',fontsize=14)
plt.ylabel('Frequency',fontsize=14)
plt.title('Histogram of N',fontsize=16)
plt.savefig(outdir + "Random_number_generator_hist_N_a_"+str(a)+".png", dpi = 300)


#plt.figure(3)
#plt.acorr(N)

plt.figure(4)
plt.acorr(U-np.mean(U))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time lag',fontsize=14)
plt.ylabel('autocorrelation function',fontsize=14)
plt.title('Autocorrelation of N',fontsize=16)
plt.savefig(outdir + "Random_number_generator_autocorr_N_a_"+str(a)+".png", dpi = 300)


#plt.figure(5)
#plt.plot(U[1:1000])