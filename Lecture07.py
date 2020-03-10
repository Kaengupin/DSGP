# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:46:51 2020

@author: r-wed
"""

import numpy as np
import matplotlib.pyplot as plt


### Settings
h, c, b  = 1, 10, 10

d = 1e-9

t  = 0.001
TU = 1
wpt = 100
T = TU/t

sigma,gamma = 1,1

dW = np.random.normal(0,1,int(T))*np.sqrt(t)


X1 = np.ones(int(T))
X2 = np.ones(int(T))
Y = np.ones(int(T))


### Functions 
def X1_val_RK(X1,X2):
    dX1 = 0
    A1, A2, A3, gamma, sigma = 0.5, 0.5, -1, 1, 1
    dX1 = (A1/gamma)*(A3*X2**2 + A2*sigma**2/(2*gamma))*X1
    return(dX1)
 
def X2_val_RK(X1,X2):
    dX2 = 0
    A1, A2, A3, gamma, sigma = 0.5, 0.5, -1, 1, 1
    dX2 = (A2/gamma)*(A3*X2**2 + A1*sigma**2/(2*gamma))*X2
    return(dX2)
    
def Y_val_RK(X1,X2,Y):
    dY = 0
    A3, gamma, epsilon = -1, 1, 0.1
    dY = A3/epsilon * X1 * X2 - gamma/epsilon**2 * Y
    return(dY)

def X1_EF(X2):
    dX1 = 0
    A1, gamma, sigma = 0.5, 1, 1
    dX1 = sigma/gamma * A1 * X2
    return(dX1)
    
def X2_EF(X1):
    dX2 = 0
    A2, gamma, sigma = 0.5, 1, 1
    dX2 = sigma/gamma * A2 * X1
    return(dX2)   
    


### Calculation
for tt in range(0,int(T)-1):
    print('Fortschritt: ' + str(int(tt*100/T)) + '%')
    ## Runge_Kutta
    i1 = X1_val_RK(X1[tt-1],X2[tt-1])
    j1 = X2_val_RK(X1[tt-1],X2[tt-1])
    k1 = Y_val_RK(X1[tt-1],X2[tt-1],Y[tt-1])
    
    i2 = X1_val_RK(X1[tt-1]+t/2*i1,X2[tt-1]+t/2*j1)
    j2 = X2_val_RK(X1[tt-1]+t/2*i1,X2[tt-1]+t/2*j1)
    k2 = Y_val_RK(X1[tt-1]+t/2*i1,X2[tt-1]+t/2*j1,Y[tt-1]+t/2*k1)
    
    i3 = X1_val_RK(X1[tt-1]+t/2*i2,X2[tt-1]+t/2*j2)
    j3 = X2_val_RK(X1[tt-1]+t/2*i2,X2[tt-1]+t/2*j2)
    k3 = Y_val_RK(X1[tt-1]+t/2*i2,X2[tt-1]+t/2*j2,Y[tt-1]+t/2*k2)
    
    i4 = X1_val_RK(X1[tt-1]+t*i3,X2[tt-1]+t*j3)
    j4 = X2_val_RK(X1[tt-1]+t*i3,X2[tt-1]+t*j3)
    k4 = Y_val_RK(X1[tt-1]+t*i3,X2[tt-1]+t*j3,Y[tt-1]+t*k3)
    
    X1[tt] = X1[tt-1] + t/6 *(i1+ 2*i2 + 2*i3 + i4) + X1_EF(X2[tt-1])*dW[tt-1]
    X2[tt] = X2[tt-1] + t/6 *(j1+ 2*j2 + 2*j3 + j4) + X2_EF(X1[tt-1])*dW[tt-1]
    Y[tt] = Y[tt-1] + t/6 *(k1+ 2*k2 + 2*k3 + k4) + (sigma/gamma)*dW[tt-1]
    

plt.figure()
plt.plot(X1)

#plt.figure()
#plt.plot(X2)

plt.figure()
plt.figure(Y)

    













