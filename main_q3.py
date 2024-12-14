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
from functionsPy.rotor import *
import pylab as plt
import os.path

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
fp = lambda x: os.path.join(inputVariables,x)

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
plt.figure()
plt.plot(wind3["t"], wind3["V_hub"])
plt.xlim(wind3["t"].min(), wind3["t"].max())
plt.xlabel("Time [s]")
plt.ylabel("Wind speeed at hub [m/s]")
plt.grid(True)

# Mean, std, max and min V_hub
print(f'Mean of V_hub is {np.mean(wind3["V_hub"])}[m/s]')
print(f'STD of V_hub is {np.std(wind3["V_hub"])}[m/s]')
print(f'Min and Max of V_hub is {np.min(wind3["V_hub"])}[m/s] & {np.max(wind3["V_hub"])}[m/s]')

# calculate Fwind from the wind time series
windForce = np.zeros_like(wind3["t"])
for i, V in enumerate(wind3["V_hub"]):
    windForce[i] += F_wind(iea22mw, wind3["V_10"], V)
wind3["F_wind"] = windForce

# plot F_wind
plt.figure()
plt.plot(wind3["t"], wind3["F_wind"]/1e6)
plt.xlabel("Time [s]")
plt.xlim(wind3["t"].min(), wind3["t"].max())
plt.ylabel("Wind forcing at hub [MN]")
plt.grid(True)
plt.show()

# calculate force statistics
print(f'Mean of F_wind is {np.mean(wind3["F_wind"])/1e6}[MN]')
print(f'STD of F_wind is {np.std(wind3["F_wind"])/1e6}[MN]')
print(f'Min and Max of F_wind is {np.min(wind3["F_wind"])/1e6}[MN] & {np.max(wind3["F_wind"])/1e6}[MN]')