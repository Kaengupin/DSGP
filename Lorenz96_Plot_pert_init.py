#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 20:23:50 2019

@author: jonpetersen
"""
import tables
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo


try:
    plt.close('all')
except:
    pass
ni=20
#Settings
path= "data/pert_init/"
title1 = "Lorenz96_XMode_RK_Pert_False_INIT_0.0"
title2 = "Lorenz96_YMode_RK_Pert_False_INIT_0.0"

#title3 = "Lorenz96_X_pert"
#title4 = "Lorenz96_Y_pert"

#fX = tables.open_file(path+title1+'.h5', mode='r')
#fY = tables.open_file(path+title2+'.h5', mode='r')
##fX1 = tables.open_file(title3+'.h5', mode='r')
##fY1 = tables.open_file(title4+'.h5', mode='r')
#
#
#X = fX.root.array_X
#Y = fY.root.array_Y


#X_cont = fX.root.array_X
#Y_cont = fY.root.array_Y

#X_pert = fX1.root.array_X
#Y_pert = fY1.root.array_Y

x=[]
xf=[]
y=[]
yf=[]
#x.append(X[:,0])
#y.append(Y[:,0])
#error = np.zeros((ni-1,2))
#error[0,0] = np.mean(abs(x[0][:]-x[0][:]))
#error[0,1] = np.mean(abs(y[0][:]-y[0][:]))

for i in range(0,ni-2):
    fX = tables.open_file(f'{path}Lorenz96_XMode_RK_Pert_True_INIT_{i}.h5', mode='r')
    fXf = tables.open_file(f'{path}Lorenz96_XMode_RK_Pert_False_INIT_{i}.h5', mode='r')
    fY = tables.open_file(f'{path}Lorenz96_YMode_RK_Pert_True_INIT_{i}.h5', mode='r')
    fYf = tables.open_file(f'{path}Lorenz96_YMode_RK_Pert_False_INIT_{i}.h5', mode='r')


    X = fX.root.array_X
    Xf = fXf.root.array_X
    Y = fY.root.array_Y
    Yf = fYf.root.array_Y

    
 
    x.append(X[:,0])
    xf.append(Xf[:,0])
    y.append(Y[:,0])
    yf.append(Yf[:,0])
    #error[i,0] = np.mean(abs(x[i][:195,:]-x[0][:195,:]))
    #error[i,1] = np.mean(abs(y[i][:195,:]-y[0][:195,:]))
    
    fX.close()
    fXf.close()
    fY.close()
    fYf.close()

#%% Pertubation Plots

X_diff_mean = []
Xdiff      = []
#
#plt.figure()
#plt.plot(x[0][:], label = "Cont")
#plt.plot(xf[0][:], label = "Pert")
#plt.legend()

for j in range(0,len(x[0])):
    for i in range(0,ni-2):
        Xdiff.append(abs(x[i][j] - xf[i][j]))
    X_diff_mean.append(np.mean(Xdiff))
    Xdiff = []

#plt.figure()
#plt.semilogy(X_diff_mean[433:700])
#
#
#X_diff = abs(X_cont[2211:2500,1]- X_pert[2211:2500,1])
#
x = np.arange(0,len(X_diff_mean[433:700]),1)
#y = 1e-9 * np.exp(0.1*x)
#
def func(v,b,c):
    return b * np.exp(c*v)

popt,pcov = spo.curve_fit(func,x[:200],X_diff_mean[433:633],bounds = ([5e-10,0.005],[2e-9,0.5]))
#
plt.figure()
plt.semilogy(x,X_diff_mean[433:700], label = "Mean Deviation")
plt.semilogy(x,func(x,popt[0],popt[1]),"r", label = r"$5\cdot 10^{-10} \cdot \mathrm{exp}(0.1126\,t)$")
#plt.title('magnitude of separation of nearby Lorenz trajectories')
plt.xlabel('time')
plt.ylabel("Value")
plt.legend(fontsize=12)
plt.savefig('Plots/Lecture_04_mean_over_initial_states.png', dpi = 300)

#





#%% Plotting

#Linienplots
#plt.figure(1)
#for pl in range(len(X[0,:])):
#    plt.plot(X[:,pl])
#
##
##plt.figure(2)
##for pl in range(len(Y[0,:])):
##    plt.plot(Y[:,pl])
#
#plt.show()
#
#
##Verteilungen
#fig=plt.figure()
#plt.hist(X[:,:].flatten(),normed=1,bins=50)
#plt.ylabel('X')
#plt.xlabel('pdf')
#plt.title('PDF of X modes')
##fig.savefig("PDF_X_Modes.pdf")    
##
##fig=plt.figure()
##plt.hist(Y[:,:].flatten(),normed=1,bins=50)
##plt.ylabel('Y')
##plt.xlabel('pdf')
##plt.title('PDF of Y modes')
##fig.savefig("PDF_Y_Modes.pdf")
##
##Hovmöller
#plt.figure()
#plt.contourf(X[:,:])
#plt.title("X-Mode")
##
##plt.figure()
##plt.contourf(Y[:,:])
##plt.title("Y-Mode")
##
##
###Punktwolke
##
##
##
###Autokorrelation
##
##
##
###Energy Cylce Terms
##E_x = []
##for i in range(len(X)):
##    E_x.append(1/2 * sum((X[i,])**2))
##
##
#
#
#fX.close()
#fY.close()
##fX1.close()
##fY1.close()



'''
x_theta = [2*np.pi /K * i for i in range(K+1)]
y_theta = [2*np.pi /(J*K) * i for i in range(K*J+1)]


plot_tt = 100
fig = plt.figure(figsize=(10, 5))
ax1 = plt.axes([0, 0, 0.4, 1], projection='polar')
#ax1 = fig.add_subplot(121, projection='polar')
ax1.plot(x_theta, X[plot_tt*10,2:-1], lw=3, zorder=10, label='X')
ax1.plot(y_theta, Y[plot_tt*10,1:-2], lw=3, label='Y')
ax1.set_rmin(-14); ax1.set_rmax(14)
l = ax1.set_rgrids([-7, 0, 7], labels=['', '', ''])[0][1]
l.set_linewidth(2)
ax1.set_thetagrids([])
ax1.set_rorigin(-22)
ax1.legend(frameon=False, loc=1);
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)



###
fig, ax = plt.subplots()
fig.set_tight_layout(True)

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))

# Plot a scatter that persists (isn't redrawn) and the initial line.
x = np.arange(0, 20, 0.1)
ax.scatter(x, x + np.random.normal(0, 3.0, len(x)))
line, = ax.plot(x, x - 5, 'r-', linewidth=2)

def update(i):
    label = 'timestep {0}'.format(i)
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    line.set_ydata(x - 5 + i)
    ax.set_xlabel(label)
    return line, ax

if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, 10), interval=200)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('line.gif', dpi=80, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()
'''
