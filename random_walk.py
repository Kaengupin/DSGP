#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:25:03 2020

@author: florian
edit: wedemannr
"""

import numpy as np
import matplotlib.pyplot as plt

try:
    plt.close('all')
except:
    pass


### Settings
dt = 2**(-6)
dt_f = 2**(-11)
n=int(1/dt_f)+1
x = np.zeros(int(1/dt_f)+1)
x_true = np.zeros(int(1/dt_f)+1)
x0 = 1.0
a = 1.5
b = 1.0
t_f = np.arange(0,1+dt_f,dt_f)
t = np.arange(0,1+dt,dt)
x_n = 1
xn = []

###Calculations

##Random Walk
dW = np.random.normal(0,1,n)*np.sqrt(dt_f)
W = np.cumsum(dW)
count = 0
##Analytical solution
for i in range(len(t_f)):
    x_true[i] = x0 * np.exp((a-0.5*b**2)*t_f[i]+b*W[i])  #Analytical 
    count += dW[i]
    if (t_f[i] in t):
        print(t_f[i])
        x_n = x_n + a*x_n*dt + b*x_n*count              #Approximation
        xn.append(x_n)
        count = 0

#%% Plot

##Random Walk
plt.figure()
plt.plot(W)
plt.title("Random Walk")

##Analytical and approximated solutions
x1 = np.linspace(0,1,len(x_true))
x2 = np.linspace(0,1,len(xn))

plt.figure()
plt.plot(x1,x_true,label="x_true")
plt.plot(x2,xn,"r",label="xn")
plt.legend()

##Plot Random Numbers dW
plt.figure()
plt.plot(dW)
plt.title("dW")
plt.show()