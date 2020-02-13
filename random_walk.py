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

dt = 2**(-4)
n=int(1/dt)+1

dW = np.random.normal(0,1,n)*np.sqrt(dt)
W = np.cumsum(dW)

plt.figure()
plt.plot(W)
plt.title("Random Walk")


x = np.zeros(int(1/dt)+1)
x_true = np.zeros(int(1/dt)+1)
x0 = 1.0
a = 1.5
b = 1.0

t = np.arange(0,1+dt,dt)

for i in range(len(t)):
    #print(i)
    x_true[i] = x0 * np.exp((a-0.5*b**2)*t[i]+b*W[i])

plt.figure()
plt.plot(x_true,label="x_true")
plt.title("x_true")



#%% Diskrete Loesung

x_n = 1
xn = []


for i in range(0,n):
    x_n = x_n + a*x_n*dt + b*x_n*dW[i]
    xn.append(x_n)

#plt.figure()
plt.plot(xn,"r",label="xn")
plt.legend()
#plt.title("xn")

plt.figure()
plt.plot(dW)
plt.title("dW")
plt.show()
