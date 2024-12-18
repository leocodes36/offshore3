'''
Filename: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_\\main_q9.py
Path: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_
Created Date: Friday, October 18th 2024, 12:13:12 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
import os
import pandas as pd
from functionsPy.common import *
import numpy as np
from functionsPy.runner import runEnvironmentalCondition
import rainflow
import matplotlib.pyplot as plt
from functionsPy.loads import *
from functionsPy.waves import *
from functionsPy.integration import *

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
TLife = 25*365*24*60*60
TSim = timeInfo['TDur'] - timeInfo["TTrans"]

# calculate the lumped sea state
scatterMat = pd.read_csv(fp("scatterMatrix.csv"))

# normalise p to be 1
scatterMat["p"] = scatterMat["p"]/np.sum(scatterMat["p"])
# calculate lifetime equivalent moment
Meq_life = (np.sum(scatterMat["p"]*scatterMat["M_eq"]**mFatigue))**(1/mFatigue) # they have different!!
print(f"lifetime equivalent moment:{Meq_life}")

"""# make waveDict
waves = loadFromJSON(fp("wave1.json"))
waves.update(timeInfo) 
waves["Hs"] = Hs
waves["Tp"] = Tp
waves["gamma"] = gamma
waves["t"] = np.arange(0., timeInfo["TDur"]+timeInfo["dt"], timeInfo["dt"])

# make windDict
wind = {}
wind["t"] = waves["t"]
wind["V_10"] = 0.
wind["V_hub"] = np.zeros_like(waves["t"])

# perform wave calculations
waves = calculateJONSWAPSpectrum(waves)
waves = generateRandomPhases(waves)
waves = calculateFreeSurfaceElevationTimeSeriesFFT(waves)
waves = calculateKinematicsFFT(waves)

plt.figure()
plt.plot(waves["eta"])

# Perform the integration
tIntegration = np.arange(0., timeInfo['TDur']+timeInfo["dt"], timeInfo["dt"])
q0 = np.array([1.0,0.])
q = ode4(dqdt, tIntegration, q0, monopile, iea22mw,
                waves, wind)

alpha = dict()
alpha["t"] = tIntegration;
alpha["alpha"] = q[:,0];
alpha["alphaDot"] = q[:,1];
alphaDotDot = np.gradient(alpha["alphaDot"], alpha["t"])
alpha["alphaDotDot"] = alphaDotDot

print("so far so good")

# calculate loads (finally)
monopile = computeElementwiseQuantities(monopile)
statLoads = calculateStaticWaveLoads(waves, monopile, alpha)
dynLoads = calculateDynamicLoads(monopile, alpha)

outputLoads = {}
outputLoads["t"] = alpha["t"]
outputLoads["F"] = statLoads["F"] + dynLoads["F"]
outputLoads["M"] = statLoads["M"] + dynLoads["M"]

# filter out the transient
Filter = outputLoads["t"] >= timeInfo["TTrans"]
outputLoads["t"] = outputLoads["t"][Filter]
outputLoads["F"] = outputLoads["F"][Filter]
outputLoads["M"] = outputLoads["M"][Filter]

plt.figure()
plt.plot(statLoads["t"], statLoads["M"], label="static")
plt.plot(dynLoads["t"], dynLoads["M"], label="dynamic")
plt.xlabel("Time [s]")
plt.ylabel("Moment [Nm]")
plt.grid()
plt.legend()

plt.figure()
plt.plot(outputLoads["t"], outputLoads["M"], label="sum")
plt.xlabel("Time [s]")
plt.ylabel("Moment [Nm]")
plt.grid()

# Rainflow count
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

# calculate equivalent moment
weighted_cycles = (amplitude**mFatigue)*(cyclesUpscaled/n_eq)
M_eq = np.sum(weighted_cycles)**(1/mFatigue)
print(M_eq)
plt.show()"""