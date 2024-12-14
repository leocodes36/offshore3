'''
Filename: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_\\main_q10.py
Path: c:\\Users\\fabpi\\OneDrive - Danmarks Tekniske Universitet\\Dokumenter\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Code\\classicalSolution\\Stud_
Created Date: Friday, October 18th 2024, 12:36:38 pm
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''

import numpy as np
import pandas as pd
from functionsPy.common import Timer, loadFromJSON
import os
import rainflow
from functionsPy.runner import *
import matplotlib.pyplot as plt

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
fp = lambda x: os.path.join(inputVariables,x)

# In this function, start from what you have in q9 and make a loop
# for all environmental conditions.

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
n_eq = 10.**7

# Time factor
TLife = 25*365*24*60*60
TSim = timeInfo['TDur'] - timeInfo["TTrans"]

# Rescale wind speed by taking into account shear factor 
shear_factor = 1/7
factor_hubheight = 2
scaleWind = factor_hubheight**shear_factor
table40 = pd.read_csv(fp("table40.csv"))
table40["V_10_scaled"] = table40["V_10"]*scaleWind

# FIX: Rescale the probabilities os that they satisfy sum(p)=1, because we left out windspeeds under 3 and above 25 earlier
table40["p"] = table40["p"]/np.sum(table40["p"])
sum_p = np.sum(table40["p"])
print("rescaled sum of probabilities",sum_p)

# Drop values outside of the cut in / cut out
dropStates, = np.where(np.logical_or(table40["V_10_scaled"]<3.0,  table40["V_10_scaled"]>25.0))
table40.drop(dropStates, inplace=True)

# Initialize the M_eq to zero
table40["M_eq"] = 0.

# Loop over ECs
for i_, ec_ in table40.iterrows():
    with Timer(f"Current EC: V={ec_['V_10']}"):
                                     
        # Here, it is more conveniente to build the dictionaries on the fly
        # rather than reading them from file.
                                               
        # Build wind dict
        wind_ = dict()
        wind_["V_10"] = ec_["V_10_scaled"]
        wind_["l"] = 340.2
        wind_["I"] = ec_["I_norm"] / 100.            
        # Make repeatable random seed, but different between runs
        wind_["randomSeed"] = i_*100 + 11

        # Build waves dict
        waves_ = dict()
        waves_["Hs"] = ec_["Hs"]
        waves_["Tp"] = ec_["Tp"]
        waves_["gamma"] = ec_["gamma_fat"]
        waves_["h"] = 34.
        waves_["z"] = np.linspace(-34.0, 0., 35)
        # Make repeatable random seed, but different between runs
        waves_["randomSeed"] = i_*100 + 22           
        
        # Run the loads calculation
        outputLoads = runEnvironmentalCondition(wind_,
                            waves_,
                            iea22mw,
                            monopile,
                            timeInfo
                            )["total"]
        
        # Then, remove transient and compute the fatigue for each sea state.                 
        # Remove transient
        Filter = outputLoads["t"] >= timeInfo["TTrans"]
        outputLoads["t"] = outputLoads["t"][Filter]
        outputLoads["F"] = outputLoads["F"][Filter]
        outputLoads["M"] = outputLoads["M"][Filter]

        # outputLoads = runEnvironmentalCondition(...)
        
        # Rainflow count
        rainflowCount = np.array(rainflow.count_cycles(outputLoads['M'])) # apply the rainflow.count_cycles routine
        amplitude = rainflowCount[:,0]/2 # transform range into amplitude
        cycles = rainflowCount[:,1] #
        cyclesUpscaled = cycles*(TLife/TSim)

        # calculate the equivalent moment here and save it in the dataframe
        weighted_cycles = (amplitude**mFatigue)*(cyclesUpscaled/n_eq)
        table40.loc[i_, "M_eq"] = np.sum(weighted_cycles)**(1/mFatigue)

# main_q11
M_eq_lifetime = (np.sum(table40["p"]*table40["M_eq"]**mFatigue))**(1/mFatigue)
print("lifetime equivalent moment",M_eq_lifetime)
# main_q12
# calculate Dlife
Dlife = np.sum(table40["M_eq"]**mFatigue*table40["p"])
# add Dstate column to table40 DataFrame
table40["Dstate"] = table40["M_eq"]**mFatigue*table40["p"]
table40["Drel"] = table40["Dstate"]/Dlife
# plot bar chart
plt.bar(table40["V_10_scaled"], table40["Drel"])
plt.xlabel("Scaled Wind Speed [m/s]")
plt.ylabel("Toxicity of sea states [-]")
plt.grid()
plt.show()

