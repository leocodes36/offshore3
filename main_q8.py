'''
Filename: c:\\Users\\fabpi\\modules\\46211assignmentsolution\\reports\\Report3\\2024\\solution\\python\\classicalSolution\\main_q8.py
Path: c:\\Users\\fabpi\\modules\\46211assignmentsolution\\reports\\Report3\\2024\\solution\\python\\classicalSolution
Created Date: Thursday, October 10th 2024, 10:27:58 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
from functionsPy.waves import *
from functionsPy.wind import *
from functionsPy.common import *
from functionsPy.monopile import forceIntegrate
import pylab as plt
import os.path

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
fp = lambda x: os.path.join(inputVariables,x)

waves = loadFromJSON(fp("wave1.json"))
timeInfo = loadFromJSON(fp("time.json"))
waves.update(timeInfo)

# Calculate the waves time vector
waves["t"] = np.arange(0., waves["TDur"], waves["dt"])

# Calculate the jonswap spectrum
waves = calculateJONSWAPSpectrum(waves)

randomSeedWaves = 1;
waves = generateRandomPhases(waves, seed=randomSeedWaves)

with Timer("slow waves"):
    waves = calculateFreeSurfaceElevationTimeSeries(waves)
    waves = calculateKinematics(waves)

# Now with FFT
with Timer("fast waves"):
    wavesFast = calculateFreeSurfaceElevationTimeSeriesFFT(waves)
    wavesFast = calculateKinematicsFFT(wavesFast)
   
# Load the wind info
wind = loadFromJSON(fp("wind4.json"))
wind.update(timeInfo)

randomWind = 11;

# Compute the wind time series
wind = calculateKaimalSpectrum(wind)
wind = generateRandomPhases(wind, seed=randomWind)

# Calculate the wind time vector
wind["t"] = np.arange(0., wind["TDur"], wind["dt"])

with Timer("Slow wind"):
    wind = calculateWindTimeSeries(wind)
with Timer("Fast wind"):
    windFast = calculateWindTimeSeriesFFT(wind)    



plt.figure()
plt.plot(waves["eta"] - wavesFast["eta"], '.')

plt.figure()
plt.plot(waves["u"][:,-1] - wavesFast["u"][:,-1], '.')

plt.figure()
plt.plot(wind["V_hub"] - windFast["V_hub"], '.')

plt.show()