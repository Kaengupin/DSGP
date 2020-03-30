'''
Hier entsteht der Ensamble Kalman Filter für das Lorenz96 Model
Die Idee ist ein Kontrolllauf XT aufzusetzen. Von diesem werden dann später nur
vereinzelnte Werte genommen, die die Observation XO darstellen. Später können
die auch noch mit einem weißen Rauschen manipuliert werden um eine Messungenauigkeit
zu simulieren. Das Modell für den Forecast XT ist das gleich Modell nur mit einem
anderen Forcing (hier 20 statt 18). Aus XT und XO wird die Analyse erstellt und
damit der neue Forecast XT erstellt. Dieser Prozess wird dann iteriert.
'''

import numpy as np
import tables
import random

'''
Setting
'''
# Lorenz96
K = 8                       # Number of X Modes
J = 32                      # K*J Number of Y Modes
True_Forcing = 18           # Forcing of the perfect model
Forecast_Forcing = 20       # Forcing of the model to simulate an error
h, c, b  = 1, 10, 10        # Params for Lorenz96

# Run Settings
t  = 0.0001                 # Timestep
TU = 10                     # How many Timeunits the Run should last
T = TU/t                    # Calc the totaö number of timesteps

# Saving and Ensamble Settings
save_timestep = 20           # Frequenz of timesteps in which data is written
analyse_cycle = 100          # Frequenz of timesteps in which a analysis cycle takes place. Must be greater then 1.
num_ensmable_members = 10    # Number of Ensamble Members
num_of_obs = K+J*K           # Number of Observations. Must stay inbetween 1 an K+J*K


def trend_X(X, Y, K = 8, J = 32, F = 18):
    '''
    Gibt den Trend für den X mode zurück.
    Übergabe ist der X und Y state und die Parameter K und J
    '''
    dX = np.zeros((K+3))
    for k in range(2,K+2):
        dX[k] = - X[k-1]*(X[k-2]-X[k+1]) - X[k] + F - h*c/b * np.sum(Y[(J*(k-1)):k*J])
    # applying cyclical condition
    dX[0] = dX[-3]
    dX[1] = dX[-2]
    dX[-1] = dX[2]
    return(dX)


def trend_Y(X, Y, K = 8, J = 32):
    '''
    Gibt den Trend für den X mode zurück.
    Übergabe ist der X und Y state und die Parameter K und J
    '''
    dY = np.zeros(((K*J)+3))
    for j in range(1,(K*J)+1):
        dY[j] = -c * b * Y[j+1] * (Y[j+2] - Y[j-1]) - c * Y[j] + (h * c)/b * X[int((j-1)/J)]
    # applying cyclical condition
    dY[0] = dY[-1]
    dY[-2] = dY[1]
    dY[-3] = dY[2]
    return(dY)

def Runge_Kutta(X, Y, dt, K = 8, J = 32, F = 18):
    '''
    Wendet Runge Kutta abhängig welcher Mode gewählt wurde.
    '''
    x1 = trend_X(X, Y, K, J, F)
    y1 = trend_Y(X, Y, K, J)

    x2 = trend_X(X+dt/2*x1, Y+dt/2*y1, K, J, F)
    y2 = trend_Y(X+dt/2*x1, Y+dt/2*y1, K, J)

    x3 = trend_X(X+dt/2*x2, Y+dt/2*y2, K, J, F)
    y3 = trend_Y(X+dt/2*x2, Y+dt/2*y2, K, J)

    x4 = trend_X(X+dt*x3, Y+dt*y3, K, J, F)
    y4 = trend_Y(X+dt*x3, Y+dt*y3, K, J)

    return(dt/6 * (x1 + 2*x2 + 2*x3 + x4), dt/6 * (y1 + 2*y2 + 2*y3 + y4))


'''
Init
'''
np.random.seed(1337)
XT = np.zeros((K+3)) + np.random.normal(size=K+3)
YT = np.zeros((J*K+3)) + np.random.normal(size=J*K+3)

'''
Analyse entspricht der Wahrheit im ersten Schritt.
Später noch mit weißem Rauschen versehen

'''
XA, YA = np.zeros((num_ensmable_members,K+3)) + XT, np.zeros((num_ensmable_members,J*K+3)) + YT


'''
Open files
'''
fXT = tables.open_file(f'data/EnKF_controll.h5', mode='w')
fXA = tables.open_file(f'data/EnKF_analysis.h5', mode='w')
atom = tables.Float64Atom()
array_XT = fXT.create_earray(fXT.root, 'array_XT', tables.Float64Atom(), shape =(0,K + K*J), title='Controll run, X and Y for Lorenz96')
array_XT.append(np.concatenate((XT[2:-1],YT[1:-2]),axis=None).reshape((1, K+J*K)))
array_XA = fXA.create_earray(fXA.root, 'array_XA', tables.Float64Atom(), shape =(0,K + K*J), title='Analysis, X and Y for Lorenz96')
array_XA.append(np.concatenate((XA[:,2:-1],YA[:,1:-2]),axis=1).reshape((1,num_ensmable_members, K+J*K)))


'''
Creating the Noice Array which is constant in time.
'''
epsi = []
noice = np.zeros((num_ensmable_members,num_of_obs,1))                   # init noice with extra dimention of transposing
for ens in range(num_ensmable_members):
    noice[ens] = (np.random.normal(0,1,num_of_obs)[np.newaxis]).T
    epsi.append(noice[ens] @ noice[ens].T)
R = np.mean(epsi, axis = 0)                                             # Mean Model Error (covariance matrix)
noice = np.squeeze(noice)                                               # get ride of the extra dimention
noice_mean = np.mean(noice, axis = 0)
# Random Set for Observations
bool_list = sorted([True] * num_of_obs + [False] * (K+J*K - num_of_obs), key=lambda k: random.random())


'''
Start des Loops
'''

print("Berechnung ueber " + str(T) + " Zeitschritte")

for tt in range(1,int(T)):

    if (tt % (T/100)) == 0:
        ''' Printing the Progress of the Run '''
        print('Fortschritt: ' + str(int(tt*100/T)) + '%')

    # True Update
    dXT, dYT = Runge_Kutta(XT, YT, t, K, J, True_Forcing)
    XT, YT = XT + dXT, YT + dYT ''' + epsi ''' ##############################

    # Analysis Update
    for ens in range(num_ensmable_members):
        dXA, dYA = Runge_Kutta(XA[ens], YA[ens], t, K, J , Forecast_Forcing)
        XA[ens], YA[ens] = XA[ens] + dXA, YA[ens] + dYA

    if (np.round(tt/save_timestep,4) % 1) == 0.:
        ''' Write data every save_timestep timesteps '''
        array_XT.append(np.concatenate((XT[2:-1],YT[1:-2]),axis=None).reshape((1, K+J*K)))
        array_XA.append(np.concatenate((XA[:,2:-1],YA[:,1:-2]),axis=1).reshape((1,num_ensmable_members, K+J*K)))

    if (np.round(tt/analyse_cycle,4) % 1) == 0.:
        ''' Analyse Update Cycle every analyse_cycle timesteps (when obs are there). '''
        # Create an Ensamble of Observations
        #obs = np.zeros((num_ensmable_members,num_of_obs))
        #for ens in range(num_ensmable_members):
        #    obs[ens] = np.concatenate((XT[2:-1],YT[1:-2]),axis=None)[bool_list] + noice[ens] # nur ein Subset auswählen

        vec = np.concatenate((XA[:,2:-1],YA[:,1:-2]),axis=1)        # State of all Ensamble Members
        ens_mean = np.mean(vec, axis = 0)
        P = np.zeros((K+J*K,K+J*K))
        for ens in range(num_ensmable_members):
            P = P + ((vec[ens]-ens_mean)[np.newaxis]).T @ ((vec[ens]-ens_mean)[np.newaxis])
        P = P / (num_ensmable_members-1)

        K = P @ H.T @ inv(H @ P @ H.T + R)

        analysis =

        XA = analysis[2:K+3] # needs a dopplecheck
        YA = analysis[K+5:-2] # needs a dopplecheck

fXT.close()
fXA.close()
