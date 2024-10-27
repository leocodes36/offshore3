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
TLife = n_eq/timeInfo['fHighCut']  #so number of equivalent cycles divided by frequency should give life
TSim = TLife * timeInfo['TDur']

# Rescale wind speed by taking into account shear factor 
# FIXME: include the correct scaleWind parameter using a shear factor of 1/7 and factor 2 on hub height
shear_factor = 1/7
factor_hubheight = 2
scaleWind = factor_hubheight**shear_factor

# Load the environmental conditions
wind = loadFromJSON(fp("wind4.json"))    
waves = loadFromJSON(fp("wave1.json"))
wind["V_10"] = wind["V_10"]*scaleWind
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

# Rainflow count
# FIXME: do the rainflow counting of fatigue here (see slides for example)
rainflowCount = np.array(rainflow.count_cycles(outputLoads['M'])) # apply the rainflow.count_cycles routine
amplitude = rainflowCount[:,0]/2 # transform range into amplitude
cycles = rainflowCount[:,1] #
cyclesUpscaled = cycles*(TLife/TSim)

plt.figure(figsize=(10, 6))
plt.hist(amplitude, bins=30, color='red', edgecolor='black')
plt.title('Histogram of Moment Amplitudes')
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.show()

# FIXME: calculate the equivalent moment here 
weighted_cycles = (amplitude**mFatigue)*(cyclesUpscaled)
M_eq = np.sum(weighted_cycles)**(1/mFatigue)
print(M_eq)