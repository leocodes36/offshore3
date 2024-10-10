'''
Filename: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classicalSolution\\Stud_\\fft_example.py
Path: c:\\Users\\fabpi\\OneDoc\\Courses\\46211_OffshoreWindEnergy\\2024\\Module3\\Lectures\\classicalSolution\\Stud_
Created Date: Tuesday, October 8th 2024, 9:43:38 am
Author: Fabio Pierella

Copyright (c) 2024 DTU Wind and Energy Systems
'''
from functionsPy.common import *
from functionsPy.waves import *
import numpy as np

# Create a wave spectrum
waves = loadFromJSON("inputVariables/wave1.json")
time = loadFromJSON("inputVariables/time.json")
waves.update(time)
waves["t"] = np.arange(0., waves["TDur"], waves["dt"])

# Calculate the spectrum
waves = calculateJONSWAPSpectrum(waves)
waves = generateRandomPhases(waves, seed=2)

# Slow solution
# the same you have inside: calculateFreeSurfaceElevationTimeSeries
with Timer("eta slow"):
    etaSlow = np.zeros_like(waves["t"])
    for i_, t_ in enumerate(waves["t"]):
        for j_, _ in enumerate(waves["f"]):
            etaSlow[i_] += waves["amplitudeSpectrum"][j_]*np.cos(2*np.pi*waves["f"][j_]*t_ + waves["randomPhases"][j_])

# Quick solution
complexSpectrum = waves["amplitudeSpectrum"]*np.exp(1j*waves["randomPhases"])
newLength = np.rint(waves["TDur"] / waves["dt"]).astype(int)
paddedSpectrum = pad2(complexSpectrum, newLength)

with Timer("fast fft"):
    etaFast = newLength*np.real(np.fft.ifft(paddedSpectrum))

import pylab as plt
plt.figure()
plt.plot(etaFast)
plt.plot(etaSlow)

plt.figure()
plt.plot(etaFast - etaSlow)