'''
Filename: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical\\main_q1.py
Path: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical
Created Date: Monday, September 30th 2024, 12:04:08 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''

from functionsPy.waves import *
from functionsPy.wind import *
from functionsPy.common import *
import pylab as plt
from main_q1 import fp

# Load the rotor info
iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = iea22mw["DRotor"]**2 * np.pi / 4

# Load the wind info
wind3 = loadFromJSON(fp("wind3.json"))
timeQ3 = loadFromJSON(fp("time.json"))
wind3.update(timeQ3)
# Calculate the time vector
wind3["t"] = np.arange(0., wind3["TDur"], wind3["dt"])

# Compute the time series
wind3 = calculateKaimalSpectrum(wind3)
wind3 = generateRandomPhases(wind3)
wind3 = calculateWindTimeSeries(wind3)

# Plot results and check if they're reasonable
plt.plot(wind3["t"], wind3["V_hub"])
plt.xlabel("Time [s]")
plt.ylabel("Wind speeed at hub [m/s]")
plt.grid(True)
plt.show()

# calculate Fwind from the wind time series