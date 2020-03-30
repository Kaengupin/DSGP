'''
Dies ist die Plotting Routine für den Ensamble Kalman Filter
'''

import tables
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.close('all')
except:
    pass

K = 8
J = 32

# Loading file
fXT = tables.open_file('data/EnKF_controll.h5', mode='r')
XT = fXT.root.array_XT

#plt.figure('first Mode')
plt.plot(XT[:,0])

plt.figure('all X modes')
for pl in range(K):
    plt.plot(XT[:,pl])

plt.figure('all Y Modes')
for pl in range(K,J*K):
    plt.plot(XT[:,pl])

plt.show()
fXT.close()
