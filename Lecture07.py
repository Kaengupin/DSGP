# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:46:51 2020

@author: r-wed
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

### Settings
h, c, b  = 1, 10, 10

d = 1e-9

t  = 0.001
TU = 1000
wpt = 100
T = TU/t

sigma,gamma = 1,1

dW = np.random.normal(0,1,int(T))*np.sqrt(t)


X1 = np.zeros(int(T))
X2 = np.zeros(int(T))
Y = np.zeros(int(T))

np.random.seed(51)
X1[0] = np.random.normal(0,np.sqrt(t),1)
np.random.seed(74)
X2[0] = np.random.normal(0,np.sqrt(t),1)
np.random.seed(11)
Y[0]  = np.random.normal(0,np.sqrt(t),1)

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
for tt in range(1,int(T)):
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
plt.plot(np.linspace(0,int(T*t),int(T)),X1)
plt.xlabel("Time (time units)",fontsize=14)
plt.ylabel(r"$X_1$",fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig(f"Plots/Lecture_07_timeseries_X1_tu{TU}.png",dpi=300)

plt.figure()
plt.plot(np.linspace(0,int(T*t),int(T)),X2)
plt.xlabel("Time (time units)",fontsize=14)
plt.ylabel(r"$X_2$",fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig(f"Plots/Lecture_07_timeseries_X2_tu{TU}.png",dpi=300)


plt.figure()
plt.plot(np.linspace(0,int(T*t),int(T)),Y)
plt.xlabel("Time (time units)",fontsize=14)
plt.ylabel(r"$Y$",fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig(f"Plots/Lecture_07_timeseries_Y_tu{TU}.png",dpi=300)

lx=int(len(X1)/5)
ly=int(len(Y)/400)

plt.figure()
ax1=plt.acorr((X1-np.mean(X1)),maxlags=lx)
plt.close()
plt.figure()
plt.plot(ax1[0][lx:]*t,ax1[1][lx:])
plt.ylim(0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.xlabel("Time lag (time units)",fontsize=14)
plt.ylabel("Autocorrelation function",fontsize=14)
plt.savefig(f"Plots/Lecture_07_autocorr_X1_tu{TU}.png",dpi=300)


plt.figure()
ax2=plt.acorr((X2-np.mean(X2)),maxlags=lx)
plt.close()
plt.figure()
plt.plot(ax2[0][lx:]*t,ax2[1][lx:])
plt.ylim(0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.xlabel("Time lag (time units)",fontsize=14)
plt.ylabel("Autocorrelation function",fontsize=14)
plt.savefig(f"Plots/Lecture_07_autocorr_X2_tu{TU}.png",dpi=300)


plt.figure()
ay=plt.acorr((Y-np.mean(Y)),maxlags=ly)
plt.close()
plt.figure()
plt.plot(ay[0][ly:]*t,ay[1][ly:])
plt.ylim(0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.xlabel("Time lag (time units)",fontsize=14)
plt.ylabel("Autocorrelation function",fontsize=14)
plt.savefig(f"Plots/Lecture_07_autocorr_Y_tu{TU}.png",dpi=300)




plt.figure()
plt.hist(X1,bins=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("bins",fontsize=14)
plt.ylabel("absolute frequency",fontsize=14)
plt.savefig(f"Plots/Lecture_07_hist_X1_tu{TU}.png",dpi=300)

    
plt.figure()
plt.hist(X2,bins=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("bins",fontsize=14)
plt.ylabel("absolute frequency",fontsize=14)
plt.savefig(f"Plots/Lecture_07_hist_X2_tu{TU}.png",dpi=300)


plt.figure()
plt.hist(Y,bins=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("bins",fontsize=14)
plt.ylabel("absolute frequency",fontsize=14)
plt.savefig(f"Plots/Lecture_07_hist_Y_tu{TU}.png",dpi=300)










