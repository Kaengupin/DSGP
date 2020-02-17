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


## Choose which Plots
Model = True
Histogram = False
Hovmoeller = False
Y_Forcing = False
Temp_Corr = False
Spat_Corr = False
Energy_cyle = False

### Save Plots
SavePlot = False
outdir = "Plot/"


K = 8
J = 32
F = 18
h, c, b  = 1, 10, 10


time = np.arange(0.0, 2000, 1)



fX = tables.open_file('data/Lorenz96_XMode_RK_Pert_False.h5', mode='r')
fY = tables.open_file('data/Lorenz96_YMode_RK_Pert_False.h5', mode='r')

X = fX.root.array_X
Y = fY.root.array_Y


if Model:
    plt.figure()
    for i in range(len(X[0,:])):
        plt.plot(X[:,i])   
    plt.xlabel("Time")
    plt.ylabel("Value")
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_XMode_Lines.png", dpi = 300)
    plt.figure()
    for i in range(len(Y[0,:])):
        plt.plot(Y[:,i]) 
    plt.xlabel("Time")
    plt.ylabel("Value")

    if SavePlot:
        plt.savefig(outdir + "Lorenz96_YMode_Lines.png", dpi = 300)
    
if Histogram:    
    ## plot distribution of X modes
    plt.figure()
    plt.hist(X[:,:K].flatten(),density=1,bins=30)
    plt.xlabel('X')
    plt.ylabel('pdf')
    plt.title('PDF of X modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_PDF_X_modes.png", dpi = 300)
    #plt.close(fig)
    
    ## plot distribution of Y modes
    plt.figure()
    plt.hist(Y[:,:J*K].flatten(),density=1,bins=30)
    plt.xlabel('Y')
    plt.ylabel('pdf')
    plt.title('PDF of Y modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_PDF_Y_modes.png", dpi = 300)
    #plt.close(fig)

if Hovmoeller:
    ##hovmoeller diagramm
    plt.figure()
    plt.contourf(X)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Hovmoeller_X_modes.png", dpi = 300)
    plt.figure()
    plt.contourf(Y)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Hovmoeller_Y_modes.png", dpi = 300)

if Y_Forcing:
## plot X mode vs Y mode forcing
    plt.figure()
    for k in range(2,K):
        plt.plot(X[:,k],h*c/b * np.sum(Y[:,(J*(k-1)+1):k*J], axis = 1) , marker ='.' , linestyle='None',Color='k')
                        
    plt.xlabel('X_k')
    plt.ylabel(r'$\frac{hc}{b}\sum_{i=J(k-1)+1}^{kJ} Y_j$')
    plt.title('X modes vs forcing of Y modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Forcing_Y_on_X_modes.png", dpi = 300)


if Temp_Corr:
    #plot 
    plt.figure()
    ac = np.zeros(100)
    for i in range(K*J):
        ts = pd.Series(Y[:,i])
        ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K*J)
    plt.plot(time[0:100], ac)
    plt.ylabel('Autokorrelation')
    plt.xlabel('time lag')
    plt.title('Temporal Autocorrelation of Y modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_TempCorr_Y_modes.png", dpi = 300)

    #plt.savefig("Autocorrelation_time_Y_modes.pdf")
    
    
    plt.figure()
    ac = np.zeros(100)
    for i in range(K):
        ts = pd.Series(X[:,i])
        ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K)
    plt.plot(time[0:100], ac)
    plt.ylabel('Autokorrelation')
    plt.xlabel('time lag')
    plt.title('Temporal Autocorrelation of X modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_TempCorr_X_modes.png", dpi = 300)


if Spat_Corr:
    ## plot spatial autocorrelation in time of Y modes
        
    plt.figure()
    ac = np.zeros(int(J*2))
    for i,t in enumerate(time):
        ts = pd.Series(Y[i,:len(time)])
        ac = ac + np.array(list(map(ts.autocorr,range(0,J*2))))/len(time)
    plt.plot(range(0,J*2), ac)
    plt.ylabel('Autocorrelation')
    plt.xlabel('time lag')
    plt.title('Spatial Autocorrelation of Y modes')
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_SpatCorr_Y_modes.png", dpi = 300)

    
        
    ## plot spatial autocorrelation in time of X modes
        
#    plt.figure()
#    ac = np.zeros(int(K/2))
#    for i,t in enumerate(time):
#        ts = pd.Series(X[i,:K])
#        ac = ac + np.array(list(map(ts.autocorr,range(0,int(K/2)))))/len(time)
#    plt.plot(range(0,int(K/2)), ac)
#    plt.ylabel('Autocorrelation')
#    plt.xlabel('time lag')
#    plt.title('Spatial Autocorrelation of X modes')
#    if SavePlot:
#        plt.savefig(outdir + "Lorenz96_SpatCorr_X_modes.png", dpi = 300)



if Energy_cyle:
    plt.figure()
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
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Energy_E_modes.png", dpi = 300)

    
    plt.figure()
    plt.plot(dtE)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Energy_dtE.png", dpi = 300)

    
    plt.figure()
    plt.plot(E)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Energy_E.png", dpi = 300)


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