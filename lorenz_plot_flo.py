#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:22:28 2020

@author: florian
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

T=int(1e5) #Anzahl Zeitschritte
dt=0.001 #Zeitintervall
s = 2000 #Anzahl gespeicherter Zeitschritte
#time = np.arange(0.0, 20001.0, 0.001)
time = np.arange(0.0, 2000, 1)
wpt = 10         # how many timesteps per time unit will be saved
tu=T*dt


fX = tables.open_file('Lorenz96_X.h5', mode='r')
fY = tables.open_file('Lorenz96_Y.h5', mode='r')

X = fX.root.array_X
Y = fY.root.array_Y

plt.figure(1)
for i in range(len(X[0,:])):
    plt.plot(X[:,i])    

plt.figure(2)
for i in range(len(Y[0,:])):
    plt.plot(Y[:,i])    
    

## plot distribution of X modes

plt.figure(3)
plt.hist(X[:,:K].flatten(),density=1,bins=30)
plt.xlabel('X')
plt.ylabel('pdf')
plt.title('PDF of X modes')
plt.savefig("pdf_X_modes.pdf")
#plt.close(fig)

## plot distribution of Y modes
    
plt.figure(4)
plt.hist(Y[:,:J*K].flatten(),density=1,bins=30)
plt.xlabel('Y')
plt.ylabel('pdf')
plt.title('PDF of Y modes')
plt.savefig("pdf_Y_modes.pdf")
#plt.close(fig)

#hovmoeller diagramm
plt.figure(5)
plt.contourf(X)

plt.figure(6)
plt.contourf(Y)


## plot X mode vs Y mode forcing
 
plt.figure(7)
for k in range(2,K):
    plt.plot(X[:,k],h*c/b * np.sum(Y[:,(J*(k-1)+1):k*J], axis = 1) , marker ='.' , linestyle='None',Color='k')
                    
plt.xlabel('X_k')
plt.ylabel(r'$\frac{hc}{b}\sum_{i=J(k-1)+1}^{kJ} Y_j$')
plt.title('X modes vs forcing of Y modes')
plt.savefig("Forcing_Y_on_X_modes.pdf")

#plot 
plt.figure(8)
ac = np.zeros(100)
for i in range(K*J):
    ts = pd.Series(Y[:,i])
    ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K*J)
plt.plot(time[0:100], ac)
plt.ylabel('Autokorrelation')
plt.xlabel('time lag')
plt.title('Temporal Autocorrelation of Y modes')
plt.savefig("Autocorrelation_time_Y_modes.pdf")


plt.figure(9)
ac = np.zeros(100)
for i in range(K):
    ts = pd.Series(X[:,i])
    ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K)
plt.plot(time[0:100], ac)
plt.ylabel('Autokorrelation')
plt.xlabel('time lag')
plt.title('Temporal Autocorrelation of X modes')
plt.savefig("Autocorrelation_time_X_modes.pdf")

  
## plot spatial autocorrelation in time of Y modes
    
#plt.figure(10)
#ac = np.zeros(int(J*2))
#for i,t in enumerate(time):
#    ts = pd.Series(Y[i,:len(time)])
#    ac = ac + np.array(list(map(ts.autocorr,range(0,J*2))))/len(time)
#plt.plot(range(0,J*2), ac)
#plt.ylabel('Autocorrelation')
#plt.xlabel('time lag')
#plt.title('Spatial Autocorrelation of Y modes')
#plt.savefig("Autocorrelation_spatial_Y_modes.pdf")

    
## plot spatial autocorrelation in time of X modes
    
plt.figure(11)
ac = np.zeros(int(K/2))
for i,t in enumerate(time):
    ts = pd.Series(X[i,:K])
    ac = ac + np.array(list(map(ts.autocorr,range(0,int(K/2)))))/len(time)
plt.plot(range(0,int(K/2)), ac)
plt.ylabel('Autocorrelation')
plt.xlabel('time lag')
plt.title('Spatial Autocorrelation of X modes')
plt.savefig("Autocorrelation_spatial_X_modes.pdf")



plt.figure(12)
E_X = []
E_Y = []
E = []
M_X = []
M_Y = []
M = []
dtE = []
for i in range(len(X[:,0])):
    ex = 1/2 * np.sum(X[i,:]**2)
    ey = 1/2 * np.sum(Y[i,:]**2)
    e = ex + ey
    mx = np.sum(X[i,:])
    my = np.sum(Y[i,:])
    m = mx + my
    E_X.append(ex)
    E_Y.append(ey)
    M_X.append(mx)
    M_Y.append(my)
    E.append(e)
    M.append(m)
    dtE.append(-ex - c*ey + F * mx)

plt.plot(E_X)
plt.plot(E_Y)
plt.plot(M_X)
plt.plot(M_Y)

plt.figure(13)
plt.plot(dtE)

plt.figure(14)
plt.plot(E)

#
## use only second half
#acor=np.zeros(len(X[0]))
#for i in range(len(X[0])):
#    norm = np.sum(X[i]**2)
#    acor[i] = np.correlate(X[:,i], X[:,i], "same")/norm
#    #acor = acor[len(acor)/2:]
#    plt.plot(acor[i])
#    #plt.acorr(X[:,i])
#    
plt.show()

fX.close()
fY.close()