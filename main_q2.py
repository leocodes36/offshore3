'''
Filename: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical\\main_q2.py
Path: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical
Created Date: Monday, September 30th 2024, 2:22:21 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''

from functionsPy.rotor import *
from functionsPy.common import *
import numpy as np
import pylab as plt
from main_q1 import fp

iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = iea22mw["DRotor"]**2 * np.pi / 4

thrust = dict()
thrust["V"] = np.arange(3., 26.)
thrust["T"] = np.zeros_like(thrust["V"])

# ASSUMPTION: V_hub is equal to V
for i_, V_ in enumerate(thrust["V"]):
    thrust["T"][i_] = F_wind(iea22mw, V_, V_)   

plt.figure()
plt.plot(thrust["V"], thrust["T"]/1e6)
plt.axvline(16, color="black", linestyle="--")
plt.xlabel("Wind speed [m/s]")
plt.ylabel("Thrust force [MN]")
plt.xlim(thrust["V"].min(), thrust["V"].max())
plt.grid(True)
plt.show()