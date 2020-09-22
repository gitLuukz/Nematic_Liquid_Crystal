import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import functions as f
import plotfunctions as pf
"""
Initialisation project
"""

epsilon = 1
T = 1
mu = 0
k = 1
delta = 1#np.pi/3.5 
#deltaN = np.flip(np.array([2,1.5,1.2,0.9,0.7,0.5,0.3]))*2
width = 24
length = 24
M = length*width #Number of molecules
N = 1000# later used for iterations of the system
NT = 25
# Trange = np.linspace(1,2,2)
# Trange = np.linspace(0.5,7,NT)

kT = np.linspace(0.1,4,NT)
# Trange = np.array([0.3,0.6,0.8,1,1.2,1.3,1.4,1.5,1.6,1.8,2.1,2.4,2.9,4,5,6,7])
# initialise locations depending on width and length
theta = np.zeros((length,width,N))
Uiold = np.zeros((length,width,N+1))
Uinew = np.zeros((length,width,N))
Ui = np.zeros((length,width,N))
UT = np.zeros((N+1,NT))
Cv = np.zeros((N,NT))
acceptance_r = np.zeros((N,NT))
T2 = np.zeros((N,NT))
T1 = np.zeros((N,NT))



# R,dx,dy = f.grid(width, length, N, theta)

# fignr = 0
# pf.orientation(fignr, R[:,:,0], dx, dy, width, length)

#%% Algorithm

theta = np.zeros((width,length,N))
[theta, UTb,avUT,scaledenergy,acceptance_r,Cv,T2,T1,Ui] = f.update(Uiold,Uinew,Ui,UT,theta,length,width,epsilon,T,mu,N,M,delta,k,acceptance_r,Cv,NT,kT,T2,T1)

#%%
thetadirec = np.linspace(-np.pi,np.pi,360)
Uidirec = -epsilon*(np.cos(2*thetadirec)+mu*np.cos(thetadirec))


#%% Fitting specific heat


x = np.zeros((int(N/2)-1,NT)) # x is confusing but it was meant as xi (and it is actually y data to make it even better :p)
taucv = np.zeros(NT)
sigmacv = np.zeros(NT)
for i in range(NT):
    x[:,i] = f.error(Cv[int(N/2):,i])
    taucv[i],pp = curve_fit(f.function,np.linspace(1,N/2-1,int(N/2-1)),x[:,i])
    sigmacv[i] = np.sqrt(2*taucv[i]/(N/2) * (np.mean(Cv[int(N/2):,i]**2) - np.mean(Cv[int(N/2):,i])**2)) 

pf.specificheat(kT,Cv[-1,:],3,sigmacv)

x = np.zeros((int(N/2)-1,NT)) # x is confusing but it was meant as xi (and it is actually y data to make it even better :p)
tauap = np.zeros(NT)
sigmaap = np.zeros(NT)
for i in range(NT):
    x[:,i] = f.error(T2[int(N/2):,i])
    tauap[i],pp = curve_fit(f.function,np.linspace(1,N/2-1,int(N/2-1)),x[:,i])
    sigmaap[i] = np.sqrt(2*tauap[i]/(N/2) * (np.mean(T2[int(N/2):,i]**2) - np.mean(T2[int(N/2):,i])**2)) 

pf.apolar(kT,T2[-1,:],6,sigmaap)


x = np.zeros((int(N/2)-1,NT)) # x is confusing but it was meant as xi (and it is actually y data to make it even better :p)
taup = np.zeros(NT)
sigmap = np.zeros(NT)
for i in range(NT):
    x[:,i] = f.error(T1[int(N/2):,i])
    taup[i],pp = curve_fit(f.function,np.linspace(1,N/2-1,int(N/2-1)),x[:,i])
    sigmap[i] = np.sqrt(2*taup[i]/(N/2) * (np.mean(T1[int(N/2):,i]**2) - np.mean(T1[int(N/2):,i])**2)) 

# pf.polar(kT,T1[-1,:],7,sigmap)

pf.both(kT, T1[-1,:], sigmap,T2[-1,:], sigmaap, 4)

x = np.zeros((int(N/2)-1,NT)) # x is confusing but it was meant as xi (and it is actually y data to make it even better :p)
tauU = np.zeros(NT)
sigmaU = np.zeros(NT)
for i in range(NT):
    x[:,i] = f.error(UTb[int(N/2+1):,i])
    tauU[i],pp = curve_fit(f.function,np.linspace(1,N/2-1,int(N/2-1)),x[:,i])
    sigmaU[i] = np.sqrt(2*tauU[i]/(N/2) * (np.mean(UTb[int(N/2):,i]**2) - np.mean(UTb[int(N/2):,i])**2)) 

# pf.polar(kT,UTb[-1,:],8,sigmaU)
pf.Energy(kT, UTb[-1,:], sigmaU, 8)

pf.energy_time(UTb[:,int(NT/2)], 9)

