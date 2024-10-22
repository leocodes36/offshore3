import numpy as np
from functionsPy.integration import ode4, dqdt
from functionsPy.common import loadFromJSON, saveToJSON
import os.path
import pylab as plt
import pandas as pd

# Location of input files
# Shorten the imports
inputVariables = "inputVariables"
savedStates = "savedStates"
fp = lambda x: os.path.join(inputVariables,x)
ss = lambda x: os.path.join(savedStates,x)

tIntegration = np.arange(0., 660., 0.1)

# Load the dictionaries
monopileDict = loadFromJSON(fp("monopile.json"))
iea22mw = loadFromJSON(fp("iea22mw.json"))
iea22mw["ARotor"] = iea22mw["DRotor"]**2*np.pi*0.25
waves5 = loadFromJSON(ss("waves4.json"))
wind5 = loadFromJSON(ss("wind4.json"))

# Perform the integration
q0 = np.array([1.0,0.])
q = ode4(dqdt, tIntegration, q0, monopileDict, iea22mw,
                waves5, wind5)

phiNodalTop = monopileDict["phiNodal"][-1]
xTTop = q[:,0]*phiNodalTop
xDotTTop = q[:,1]*phiNodalTop
f,ax = plt.subplots(2, sharex=True)
ax[0].plot(tIntegration, xTTop)
ax[1].plot(tIntegration, xDotTTop)
ax[0].grid();ax[1].grid()
ax[1].set_xlabel("t[s]")
ax[0].set_ylabel("XTTop[m]")
ax[1].set_ylabel("XTDotTop[m/s]")
plt.show()

alpha = dict()
alpha["t"] = tIntegration;
alpha["alpha"] = q[:,0];
alpha["alphaDot"] = q[:,1];

plt.plot(alpha["t"], alpha["alpha"])
plt.xlabel("Time [s]")
plt.ylabel("alpha_1")
plt.grid()
plt.show()

DataFrame = pd.DataFrame.from_dict(alpha)
print(DataFrame.describe())

saveToJSON(alpha, ss("q.json"))