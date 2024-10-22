'''
Filename: c:\\Users\\fabpi\\modules\\46211assignmentsolution\\reports\\Report3\\2024\\solution\\python\\classicalSolution\\main_q7.py
Path: c:\\Users\\fabpi\\modules\\46211assignmentsolution\\reports\\Report3\\2024\\solution\\python\\classicalSolution
Created Date: Monday, October 7th 2024, 11:52:58 am
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
import pylab as plt
from functionsPy.common import loadConstants, loadFromJSON, downsample, saveToJSON    
import os.path
import numpy as np
from functionsPy.loads import *
from functionsPy.monopile import *

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
savedStates = "savedStates"
fp = lambda x: os.path.join(inputVariables,x)
ss = lambda x: os.path.join(savedStates,x)

# Load the structural motions
q = loadFromJSON(ss("q.json"))
alphaDotDot = np.gradient(q["alphaDot"], q["t"])
q["alphaDotDot"] = alphaDotDot

# Load the wind and waves
wind = loadFromJSON(ss("wind4.json"))
waves = loadFromJSON(ss("waves4.json"))

# Load the monopile and the rotor
iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = iea22mw["DRotor"]**2*np.pi/4

monopileDict = loadFromJSON(fp("monopile.json"))

# Calculate the static wind loads
windDownsampled = downsample(wind, dropEvery=2, listOfFields=["t", "V_hub"])
windLoads = calculateStaticWindLoads(windDownsampled, iea22mw, 
                                monopileDict, q)

# Calculate the static wave loads
wavesDownsampled = downsample(waves, dropEvery=2, listOfFields=["t", "u", "ut"])
waveLoads = calculateStaticWaveLoads(wavesDownsampled, monopileDict, q)

# Calculate the dynamic loads
monopileDict = computeElementwiseQuantities(monopileDict)
dynamicLoads = calculateDynamicLoads(monopileDict, q)

saveToJSON(dynamicLoads, "savedStates/dynLoads.json")

# Figure
plt.figure()
plt.plot(dynamicLoads["t"], dynamicLoads["M"], label="MDyn")
plt.plot(windLoads["t"], windLoads["M"], label="Wind")
plt.plot(waveLoads["t"], waveLoads["M"], label="Waves")
plt.legend()
plt.grid()
plt.show()