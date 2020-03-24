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
Model = False
Histogram = True
Histogram_Diff = True
Hovmoeller = False
Y_Forcing = False
Temp_Corr = False
Spat_Corr = False
Energy_cyle = False

### Save Plots
SavePlot = True
outdir = "Plots/"

t=0.0001

K = 8
J = 32
F = 18
h, c, b  = 1, 10, 10
TU = 100


time = np.arange(0.0, 2000, 1)


#fX = tables.open_file('data/Lorenz96_XMode_RRK_Pert_False.h5', mode='r')
#fY = tables.open_file('data/Lorenz96_YMode_RRK_Pert_False.h5', mode='r')

fX = tables.open_file('data/Lorenz96_XMode_RK_Pert_False_TU_100.h5', mode='r')
fY = tables.open_file('data/Lorenz96_YMode_RK_Pert_False_TU_100.h5', mode='r')

fXr = tables.open_file('data/Lorenz96_XMode_RRK_Pert_False_TU_100.h5', mode='r')
fYr = tables.open_file('data/Lorenz96_YMode_RRK_Pert_False_TU_100.h5', mode='r')


X = fX.root.array_X
Y = fY.root.array_Y

Xr = fXr.root.array_X
Yr = fYr.root.array_Y

#%%

if Model:
    plt.figure()
    for i in range(len(X[0,:])):
        plt.plot(np.linspace(0,TU,len(X[:])),X[:,i])   
    plt.xlabel("Time Unit")
    plt.ylabel("Value")
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_XMode_Lines_TU_"+str(TU)+".png", dpi = 300)
    plt.figure()
    for i in range(len(Y[0,:])):
        plt.plot(np.linspace(0,TU,len(X[:])),Y[:,i]) 
    plt.xlabel("Time Unit")
    plt.ylabel("Value")

    if SavePlot:
        plt.savefig(outdir + "Lorenz96_YMode_Lines_TU_"+str(TU)+".png", dpi = 300)
    
if Histogram:    
    ## plot distribution of X modes
    plt.figure()
    plt.hist(X[:,:K].flatten(),bins=36,density=1,range=(-15,20))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('X',fontsize=14)
    plt.ylabel('relative frequency',fontsize=14)
    plt.title('PDF of X modes',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_PDF_X_modes_TU_"+str(TU)+".png", dpi = 300)
    #plt.close(fig)
    
    plt.figure()
    plt.hist(Xr[:,:K].flatten(),bins=36,density=1,range=(-15,20))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('X',fontsize=14)
    plt.ylabel('relative frequency',fontsize=14)
    plt.title('PDF of X modes with Wilks scheme',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_PDF_Xr_modes_TU_"+str(TU)+".png", dpi = 300)
    #plt.close(fig)
    
    ## plot distribution of Y modes
    #plt.figure()
    #plt.hist(Y[:,:J*K].flatten(),density=1,bins=30)
    #plt.xticks(fontsize=12)
    #plt.yticks(fontsize=12)
    #plt.xlabel('Y',fontsize=14)
    #plt.ylabel('PDF',fontsize=14)
    #plt.title('PDF of Y modes',fontsize=16)
    #if SavePlot:
    #    plt.savefig(outdir + "Lorenz96_PDF_Y_modes_TU_"+str(TU)+".png", dpi = 300)
    #plt.close(fig)

if Histogram_Diff:    
    plt.figure()
    #print(Xr[:,:K].shape)
    #print(X[:,:K].shape)
    hist=plt.hist((X[:,:K].flatten()),bins=36,range=(-15,20))
    histr=plt.hist((Xr[:,:K].flatten()),bins=36,range=(-15,20))
    plt.close()
    plt.figure()
    hist_diff=histr[0]-hist[0]
    plt.bar(np.linspace(-15,20,36),hist_diff/np.sum(hist[0]),width=1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(-15,20)
    plt.xlabel('X',fontsize=14)
    plt.ylabel('difference in relative frequency',fontsize=14)
    plt.title('Difference to X modes in deterministic scheme',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_PDF_X_modes_diff_TU_"+str(TU)+".png", dpi = 300)

if Hovmoeller:
    ##hovmoeller diagramm
    plt.figure()
    plt.contourf(X[9000:10000])
    print(X.shape)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("X-Modes",fontsize=14)
    plt.ylabel("Time",fontsize=-14)
    plt.title('Hovmoeller diagram of X modes',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Hovmoeller_X_modes_TU_"+str(TU)+".png", dpi = 300)
    plt.figure()
    plt.contourf(Y[9000:10000,:50])
    print(Y[9900:10000,:50].shape)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Y-Modes",fontsize=14)
    plt.ylabel("Time",fontsize=14)
    plt.title('Hovmoeller diagram of Y modes',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Hovmoeller_Y_modes_TU_"+str(TU)+".png", dpi = 300)

if Y_Forcing:
## plot X mode vs Y mode forcing
    plt.figure()
    for k in range(2,K):
        plt.plot(X[:,k],h*c/b * np.sum(Y[:,(J*(k-1)+1):k*J], axis = 1) , marker ='.' , linestyle='None',Color='k')
                        
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('X_k',fontsize=14)
    plt.ylabel(r'$\frac{hc}{b}\sum_{i=J(k-1)+1}^{kJ} Y_j$',fontsize=14)
    plt.title('X modes vs forcing of Y modes',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Forcing_Y_on_X_modes_TU_"+str(TU)+".png", dpi = 300)


if Temp_Corr:
    #plot 
    plt.figure()
    ac = np.zeros(100)
    for i in range(K*J):
        ts = pd.Series(Y[:,i])
        ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K*J)
    plt.plot(time[0:100], ac)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Autokorrelation',fontsize=14)
    plt.xlabel('time lag',fontsize=14)
    plt.title('Temporal Autocorrelation of Y modes',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_TempCorr_Y_modes_TU_"+str(TU)+".png", dpi = 300)

    #plt.savefig("Autocorrelation_time_Y_modes.pdf")
    
    
    plt.figure()
    ac = np.zeros(100)
    for i in range(K):
        ts = pd.Series(X[:,i])
        ac = ac + np.array(list(map(ts.autocorr,range(0,100))))/(K)
    plt.plot(time[0:100], ac)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Autokorrelation',fontsize=14)
    plt.xlabel('time lag',fontsize=14)
    plt.title('Temporal Autocorrelation of X modes with Wilks scheme',fontsize=16)
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_TempCorr_X_modes_TU_"+str(TU)+".png", dpi = 300)


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
        plt.savefig(outdir + "Lorenz96_SpatCorr_Y_modes_TU_"+str(TU)+".png", dpi = 300)

    
        
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
    
    plt.figure()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.plot(E_X,label='E_X')
    plt.plot(E_Y,label='E_Y')
    plt.title('Energy terms',fontsize=16)
    plt.legend()
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Energy_E_modes_TU_"+str(TU)+".png", dpi = 300)


    plt.figure()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.plot(M_X,label='M_X')
    plt.plot(M_Y,label='M_Y')
    plt.title('Momentum terms',fontsize=16)
    plt.legend()
    if SavePlot:
        plt.savefig(outdir + "Lorenz96_Momentum_M_modes_TU_"+str(TU)+".png", dpi = 300)

    
    plt.figure()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.plot(dtE)
    #if SavePlot:
    #    plt.savefig(outdir + "Lorenz96_Energy_dtE.png_TU_"+str(TU)+".png", dpi = 300)

    
    plt.figure()
    plt.plot(E)
    #if SavePlot:
    #    plt.savefig(outdir + "Lorenz96_Energy_E.png_TU_"+str(TU)+".png", dpi = 300)


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

fXr.close()
fYr.close()