'''
Filename: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classicalSolution\\Stud_\\ode_example.py
Path: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classicalSolution\\Stud_
Created Date: Tuesday, October 8th 2024, 9:20:30 am
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
from numpy import array, arange, pi, sin 
from functionsPy.integration import ode4
import matplotlib.pyplot as plt

# This is an example implementation.
# Initial conditions
y0 = array([0,1])

# Generalized masses and stiffnesses
GM = 1.0
GD = 0.1
GK = 1.0

# Generalized forcing
omega_f = 2*pi
A_f = 0.0
phi_f = 0.

GF = lambda t: A_f*sin(omega_f*t+phi_f)

dqdt = lambda t,q: array([
        q[1],
        (-GK*q[0] -GD*q[1] + GF(t))/GM])

tspan = arange(0., 10., 0.1)

# Integration 
y = ode4(dqdt, tspan, y0)

# Plotting
plt.figure()
plt.plot(tspan, y[:,0], label = 'position')
plt.plot(tspan, y[:,1], label = 'velocity')
plt.legend()
plt.grid()
plt.xlabel('time[s]')
plt.show()
