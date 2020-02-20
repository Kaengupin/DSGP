#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:50:59 2020

@author: u301023
"""

import numpy as np
import matplotlib.pyplot as plt
import tables
import pandas as pd

try:
    plt.close('all')
except:
    pass

F=20
h=1
c=1
b=1
K=8
J=6

nt  = 10                                    # number of time steps
dt  = [1e-4,2e-4,4e-4,8e-4,1.6e-3,3.2e-3,6.4e-3,1.28e-2,2.56e-2]   # vector of different time steps
TU  = 1                                     # Time unit
wpt = 100                                   # how many timesteps per time unit will be saved

title1 = "Lorenz96_XMode"
title2 = "Lorenz96_YMode"


### choose scheme
scheme = "RK"     # EF:Euler Forward, RK:Runge Kutta, RRK:Reduced Runge Kutta

# fX = tables.open_file(f'data/TS/{title1}_{scheme}_TS_1.h5', mode='r')
# fY = tables.open_file(f'data/TS/{title2}_{scheme}_TS_1.h5', mode='r')

# X = fX.root.array_X
# Y = fY.root.array_Y

# x0 = X[:]
# y0 = Y[:]

# fX.close()
# fY.close()

error = np.zeros((nt-1,2))

x = []
y = []

for i in range(0,nt-1):
    fX = tables.open_file(f'data/TS/{title1}_{scheme}_TS_{i+1}.h5', mode='r')
    fY = tables.open_file(f'data/TS/{title2}_{scheme}_TS_{i+1}.h5', mode='r')

    X = fX.root.array_X
    Y = fY.root.array_Y
 
    x.append(X[:])
    y.append(Y[:])
    error[i,0] = np.mean(abs(x[i][:195,:]-x[0][:195,:]))
    error[i,1] = np.mean(abs(y[i][:195,:]-y[0][:195,:]))
    
    fX.close()
    fY.close()

plt.figure()
plt.loglog(dt,error[:,0])

plt.figure()
plt.loglog(dt,error[:,1])