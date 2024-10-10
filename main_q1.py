'''
Filename: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical\\main_q1.py
Path: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classical
Created Date: Monday, September 30th 2024, 12:04:08 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''

from functionsPy.waves import *
from functionsPy.common import *
from functionsPy.monopile import forceIntegrate
import pylab as plt
import os.path

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
fp = lambda x: os.path.join(inputVariables,x)

# Question 1
wavesQ1 = loadFromJSON(fp("wave1.json"))
    
# Load the time discretization info
timeQ1 = loadFromJSON(fp("time.json"))
wavesQ1.update(timeQ1)

# Calculate the time vector
wavesQ1["t"] = np.arange(0., wavesQ1["TDur"], wavesQ1["dt"])

# Calculate the jonswap spectrum
wavesQ1 = calculateJONSWAPSpectrum(wavesQ1)

randomSeedWaves = 1;
wavesQ1 = generateRandomPhases(wavesQ1, seed=randomSeedWaves)

wavesQ1 = calculateFreeSurfaceElevationTimeSeries(wavesQ1)
wavesQ1 = calculateKinematics(wavesQ1)

monopileDict = loadFromJSON(fp("monopile.json"))

forceQ1 = dict()
forceQ1["t"] = wavesQ1["t"]
forceQ1["F"], forceQ1["M"] = np.zeros_like(wavesQ1["t"]), np.zeros_like(wavesQ1["t"])

for i_, t_ in enumerate(wavesQ1["t"]):
    forceQ1["F"][i_], forceQ1["M"][i_]  = forceIntegrate(monopileDict, wavesQ1["u"][i_,:], wavesQ1["ut"][i_,:],
        wavesQ1["z"], 0.)

# check spectrum is correct: for 4*sigma_eta vs Hs
print("4 sigma_eta:", 4*np.std(wavesQ1["eta"]))
print("Hs", wavesQ1["Hs"])

# plot result for eta
plt.figure()
plt.plot(wavesQ1["t"], wavesQ1["eta"])
plt.xlabel("Time [s]")
plt.ylabel("Free surface elevantion [m]")
plt.xlim(0, wavesQ1["TDur"])
plt.grid(True)

# plot result for force
plt.figure()
plt.plot(forceQ1["t"], forceQ1["F"]/1e6)
plt.xlabel("Time [s]")
plt.ylabel("Force [MN]")
plt.xlim(0, wavesQ1["TDur"])
plt.grid(True)

# plot result for force
plt.figure()
plt.plot(forceQ1["t"], forceQ1["M"]/1e6)
plt.xlabel("Time [s]")
plt.ylabel("Moment [MN]")
plt.xlim(0, wavesQ1["TDur"])
plt.grid(True)

# show figures
plt.show()
