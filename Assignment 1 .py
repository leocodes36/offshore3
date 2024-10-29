# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:22:54 2024

@author: markj
"""
import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as sopt
from functions import forceIntegrate

## input parameters and monopileData as "structure"
h = 34. #height
g = 9.81 #gravity
dt = 0.1 #time steps
Nz = 40 #number of arrays

monopileData = dict()

monopileData["CM"] = 2. #center of mass
monopileData["CD"] = 1. #drag coefficient
monopileData["D"] = 10. #diameter
monopileData["h"] = h #height of monopile
monopileData["rho"] = 1025. #density

if True:
    
    # regular waves
    H=np.array([6]);
    T=np.array([12]); #period
    f=1/T;
    a=H/2; #amplitude
    epsilon=np.array([0]); #free surface elevation
    TDur=2*T[0]; #duration
    
else:
    # irregular waves
    
    Hs = 11.3
    Tp = 18.5
    
    fHighCut = 0.50
    TDur = 600.
    df = 1/TDur
    f = np.arange(df, fHighCut+df, df)
    
    # code the JONSWAP spectrum here. Right now it is just a uniform
    # spectrum.
    S=np.ones(np.size(f));
    
    # Plot spectrum
    plt.figure()
    plt.plot(f, S)
    plt.xlabel(r"$f[Hz]$")
    plt.ylabel(r"$S_\eta[m^2/Hz]$")
    
    # Calc a and assign random phases
    np.random.seed(22)
    a = np.sqrt(2*S*df)
    epsilon = np.random.rand(len(a))*2*np.pi
    
## Solve the dispersion relation to get the wave numbers k

N = len(f)
omega = 2*np.pi*f
k = np.zeros_like(f)  #generates array of zeros

myFun = lambda k, omega, g, h: omega**2 - g*k*np.tanh(k*h) # % modify dispersion relation to work for all depths
myFunPrime = lambda k, omega, g, h: -g*(k*h*(1-np.tanh(k*h)**2)+np.tanh(k*h));

kGuess = omega[0] / np.sqrt(g*h)

for j in range(N):
    k[j] = sopt.root_scalar(lambda x: myFun(x,omega[j],g,h), fprime=lambda x: myFunPrime(x,omega[j], g,h), x0=kGuess, method='newton').root;
    kGuess = k[j]
    
L = 2*np.pi/k    
print(L)
# L = 2*pi/k  %Q1: L=184.5 m

# Calculate free surface elevation and force in a loop
t = np.arange(0., TDur, dt)
Nt = len(t)

Eta=np.zeros(Nt)
U = np.zeros((Nz, Nt))
Ut = np.zeros((Nz, Nt))

Fx = np.zeros(Nt)
My = np.zeros(Nt)

for it in range(0,Nt):
    Eta[it] = sum(   a*np.cos(omega*t[it]+epsilon)    ) ;
    zPhys = np.linspace(-h,Eta[it],Nz);
    #zCalc = np.linspace(-h, 0., Nz) # Uuse zCalc in the calculation of velocity and acceleration to apply Wheeler stretching

    for iz in range(Nz):  #change zPhys to zCalc to apply Wheeler stretching
        U[iz, it]= np.sum(    a*omega*np.cosh(k*(zPhys[iz]+h))/np.sinh(k*h)*np.cos( omega*t[it]+epsilon) )   #calculate u[m/s] no kx as we chose x as 0
        Ut[iz, it] = np.sum ( -a*omega*np.cosh(k*(zPhys[iz]+h))/np.sinh(k*h)*np.sin(omega*t[it]+epsilon))   # differentiate with respect to x 
        
    
        

    [F, M]= forceIntegrate(monopileData, U[:,it], Ut[:,it], zPhys)
    
    Fx[it] = F
  

# Make some plots

f2,ax2 = plt.subplots(3, sharex=True)

ax2[0].plot(t,Eta,'b')
ax2[0].set_xlabel('t [s]')
ax2[0].set_ylabel('η[m]')
ax2[0].grid()

ax2[1].plot(t,Fx,'k')
ax2[1].set_xlabel('t [s]')
ax2[1].set_ylabel('Force [N]')
ax2[1].grid()

ax2[2].plot(t,My,'r')
ax2[2].set_xlabel('t [s]')
ax2[2].set_ylabel('Moment [Nm]')
ax2[2].grid()

f2.tight_layout()

# Make some statistics
print("{:8} {:^15} {:^15} {:^15}".format('', 'eta', 'Fx', 'My'))
columnFormat = "{:8} {:^15.2f} {:^15.2f} {:^15.2f}"
print(columnFormat.format('mean', np.mean(Eta), np.mean(Fx), np.mean(My)))
print(columnFormat.format('std', np.std(Eta), np.std(Fx), np.std(My)))
print(columnFormat.format('min', np.min(Eta), np.min(Fx), np.min(My)))
print(columnFormat.format('max', np.max(Eta), np.max(Fx), np.max(My)))

f3,ax3 = plt.subplots(3, sharex=True)

ax3[0].plot(t,Eta,'b', label="eta")
ax3[0].set_xlabel('t [s]')
ax3[0].set_ylabel('η[m]')
ax3[0].grid()

ax3[1].plot(t,U[0,:],'b', label="velocity u, z=-h")
ax3[1].plot(t,Eta-U[0,:],'r', label="velocity u, z=0")  # add Eta - U
ax3[1].set_xlabel(r'$t [s]$')
ax3[1].set_ylabel(r'$u [m/s]$')
ax3[1].grid()
ax3[1].legend()

ax3[2].plot(t,Ut[0,:],'b', label="acc u_t, z=-h")
ax3[2].plot(t,Eta-Ut[0,:],'r', label="acc u_t, z=0")    # Add Eta -
ax3[2].set_xlabel(r'$t [s]$')
ax3[2].set_ylabel(r'$u_t [m/s^2]$')
ax3[2].grid()
ax3[2].legend()

f3.tight_layout()
