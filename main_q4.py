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
from functionsPy.monopile import *
from functionsPy.rotor import *
import pylab as plt
import os.path

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
savedStates = "savedStates"
fp = lambda x: os.path.join(inputVariables,x)
ss = lambda x: os.path.join(savedStates,x)

# Load the wind info
wind4 = loadFromJSON(fp("wind4.json"))
time4 = loadFromJSON(fp("time.json"))
wind4.update(time4)

# Load the waves stuff
waves4 = loadFromJSON(fp("wave1.json"))
waves4.update(time4)

wind4["t"] = np.arange(0., wind4["TDur"], wind4["dt"])
waves4["t"] = np.arange(0., waves4["TDur"], waves4["dt"])

# Compute the time series
wind4 = calculateKaimalSpectrum(wind4)
wind4 = generateRandomPhases(wind4, seed=1)
wind4 = calculateWindTimeSeries(wind4)

# Compute waves
waves4 = calculateJONSWAPSpectrum(waves4)
waves4 = generateRandomPhases(waves4, seed=2)
waves4 = calculateFreeSurfaceElevationTimeSeries(waves4)
waves4 = calculateKinematics(waves4)
h = waves4["h"]

# Load the monopile
monopileDict = loadFromJSON(fp("monopile.json"))

# Wave force
waveForce = dict()
waveForce["t"] = waves4["t"]
waveForce["F"], waveForce["M"] = np.zeros_like(waves4["t"]), np.zeros_like(waves4["t"])

for i_, t_ in enumerate(waves4["t"]):
	# FIXME: call the forceIntegrate function to get the wave loads
    waveForce["F"][i_], waveForce["M"][i_]  = 0., 0.
    
# Wind force

iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = iea22mw["DRotor"]**2*np.pi/4

windForce = dict()
windForce["t"] = wind4["t"]
windForce["F"], windForce["M"] = np.zeros_like(wind4["t"]), np.zeros_like(wind4["t"])

for i_, t_ in enumerate(wind4["t"]):
	# FIXME: call the F_wind to get the total wind force
    windForce["F"][i_] = 0.
    windForce["M"][i_] = windForce["F"][i_]*(monopileDict["zBeamNodal"][-1]+h)
    
plt.figure()
plt.plot(waves4["t"], waves4["eta"])

plt.figure()
plt.plot(waveForce["t"], waveForce["F"])

# Save for later use
os.makedirs(savedStates, exist_ok=True)
saveToJSON(waves4, ss("waves4.json"))
saveToJSON(wind4, ss("wind4.json" ))