'''
Filename: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_\\main_q9.py
Path: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_
Created Date: Friday, October 18th 2024, 12:13:12 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
import os
import pandas as pd
from functionsPy.common import loadFromJSON, Timer
import numpy as np
from functionsPy.runner import runEnvironmentalCondition
import rainflow
import matplotlib.pyplot as plt
from functionsPy.waves import *
from functionsPy.wind import *

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
fp = lambda x: os.path.join(inputVariables,x)

# Load the rotor and the monopile which do not change
iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = 0.25*np.pi*iea22mw["DRotor"]**2
iea22mw["VCutIn"] = 3.0; iea22mw["VCutOut"] = 25.0
monopile = loadFromJSON(fp("monopile.json"))

# Time information
timeInfo = dict()
timeInfo["TDur"] = 3660.
timeInfo["TTrans"] = 60.
timeInfo["dt"] = 0.05
timeInfo["fHighCut"] = 0.5

# Fatigue parameters
mFatigue = 4.
n_eq = 10**7

# Time factor
# FIXME: Calculate TLife and Tsim
# FIXED
TLife = 25*365*24*60*60
TSim = timeInfo['TDur'] - timeInfo["TTrans"]

# no rescaling, since wind4.json contains already scaled windspeed 17.66 m/s

# Load the environmental conditions
wind = loadFromJSON(fp("wind4.json"))    
waves = loadFromJSON(fp("wave1.json"))
wind["V_10"] = wind["V_10"]
wind["randomSeed"], waves["randomSeed"] = 1, 2

# Run the loads calculation
outputLoads = runEnvironmentalCondition(wind,
                    waves,
                    iea22mw,
                    monopile,
                    timeInfo
                    )["total"]
                        
# Remove transient
Filter = outputLoads["t"] >= timeInfo["TTrans"]
outputLoads["t"] = outputLoads["t"][Filter]
outputLoads["F"] = outputLoads["F"][Filter]
outputLoads["M"] = outputLoads["M"][Filter]

plt.plot(outputLoads["t"], outputLoads["M"]/1e6)
plt.xlabel("Time [s]")
plt.ylabel("Moment [MNm]")
plt.grid()

# Rainflow count
# FIXME: do the rainflow counting of fatigue here (see slides for example)
# FIXED
rainflowCount = np.array(rainflow.count_cycles(outputLoads['M'])) # apply the rainflow.count_cycles routine
amplitude = rainflowCount[:,0]/2 # transform range into amplitude
cycles = rainflowCount[:,1] #
cyclesUpscaled = cycles*(TLife/TSim)

plt.figure()
plt.hist(amplitude/1e6, bins=30, color='red', edgecolor='black')
plt.title('Histogram of Moment Amplitudes')
plt.xlabel('Amplitude [MNm]')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# plot wind and wave Spectra
wind = calculateKaimalSpectrum(wind)
waves = calculateJONSWAPSpectrum(waves)

# calculate natural frequency
omega0 = np.sqrt(monopile["GK"]/monopile["GM"])
print(f"Nat. Frequency:{omega0}")

plt.figure()
plt.plot(wind["f"], wind["Spectrum"]/np.max(wind["Spectrum"]), label="Kaimal Spectrum (Wind)")
plt.plot(waves["f"], waves["Spectrum"]/np.max(waves["Spectrum"]), label="JONSWAP Spectrum (Waves)")
plt.axvline(x=0.16, color="red", linestyle="--", linewidth=0.5, label="Natural Frequency")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Normalised Spectral Density [-]")
plt.grid()
plt.legend()

# FIXME: calculate the equivalent moment here 
# FIXED
weighted_cycles = (amplitude**mFatigue)*(cyclesUpscaled/n_eq)
M_eq = np.sum(weighted_cycles)**(1/mFatigue)
print(f"M_eq:{M_eq}")
plt.show()
