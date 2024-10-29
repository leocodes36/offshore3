import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as sopt
from functions import forceIntegrate
import common as cm

# ----------------------------------------------- Part 1 - Exercise 2.1 -----------------------------------------------
## input parameters and monopileData as "structure"
h = 34.  # depth
g = 9.81  # gravity
dt = 0.05  # time steps
Nz = 20  # number of columns

monopileData = dict()  # base structure

monopileData["CM"] = 2.  # inertia coefficient
monopileData["CD"] = 1.  # drag coefficient
monopileData["D"] = 10.  # diameter of the monopile
monopileData["h"] = h  # depth
monopileData["rho"] = 1025.  # density of water


# irregular waves
Hs = 2.19  # significant wave height
Tp = 6.37

fHighCut = 0.5
TDur = 660.
df = 1 / TDur
f = np.arange(df, fHighCut, df) # this is different in assignment 3
fp = 1/Tp
gamma = 1.0
"""
S = np.ones(np.size(f));


if f.any() <= fp:
    sigma = 0.07
else:
    sigma = 0.09
"""

# Spectral width parameter
sigma = np.ones(len(f))
sigma[f>fp] = 0.09
sigma[f<=fp] = 0.07
#Rewrite JS spectrum from Assignment 1 since wrong
#S = np.ones(np.size(f));
#S = np.ones_like(f)
S = 5/16*(Hs**2)*Tp*((f/fp)**(-5))*np.exp((-5/4)*((f/fp)**(-4)))*(1-0.287*np.log(gamma)*gamma)**np.exp((-0.5)*(((f/fp)-1)/sigma)**2)

# Plot spectrum
plt.figure()
plt.plot(f, S)
plt.xlabel(r"$f[Hz]$")
plt.ylabel(r"$S_\eta[m^2/Hz]$")
plt.grid()



# Calc a and assign random phases
# Initialize the waves dictionary
waves = {
    "Spectrum": S,  # Ensure 'Spectrum' key is present
}

# Define the LCG function to generate random phases
def lcg(seed, a=1103515245, c=12345, m=2**31, n=1):
    numbers = np.ones(n) * seed
    for i in range(1, n):
        numbers[i] = (a * numbers[i - 1] + c) % m
    return 2 * np.pi * np.array([x / m for x in numbers])  # Normalize to [0, 1]

# Function to generate random phases and include them in the waves dictionary
def generateRandomPhases(inputDict, seed=2):
    # Generate random phases
    phi = lcg(seed, n=len(inputDict["Spectrum"]))
    
    outputDict = dict()
    outputDict.update(inputDict)    
    outputDict["randomPhases"] = phi  # Add the random phases
    outputDict["epsilon"] = phi  # If epsilon is meant to be random phases, add it too

    return outputDict

# Use the function to generate random phases
randomSeedWaves = 1
waves = generateRandomPhases(waves, seed=randomSeedWaves)
a = np.sqrt(2 * S * df)

"""
#np.random.seed(2)
seed = 1
rng = np.random.default_rng(seed)
a = np.sqrt(2 * S * df)
epsilon = rng.random(len(a)) * 2 * np.pi  # Generate random phases in range [0, 2π]

#epsilon = np.random.rand(len(a)) * 2 * np.pi
"""

## Solve the dispersion relation to get the wave numbers k

N = len(f)  # length of the frequency array
omega = 2 * np.pi * f  # radian frequency
k = np.zeros_like(f)  # wave number, zero array with the same size as the frequency

# We isolate the equation to make it equal to zero
myFun = lambda k, omega, g, h: omega ** 2 - (g * k * np.tanh(k * h))
myFunPrime = lambda k, omega, g, h: -g * (
            k * h * (1 - np.tanh(k * h) ** 2) + np.tanh(k * h))  # full dispersion relation

kGuess = omega[0] / np.sqrt(g * h)  # first estimate of k

for j in range(N):
    k[j] = sopt.root_scalar(lambda x: myFun(x, omega[j], g, h), fprime=lambda x: myFunPrime(x, omega[j], g, h),
                            x0=kGuess, method='newton').root;
    kGuess = k[j]

L = 2 * np.pi / k  # %Q1: L=184.5 m
#print("Part 1 - Exercise 1 = " + str(L))

# Calculate free surface elevation and force in a loop
t = np.arange(0., TDur,
              dt)  # generates an array of evenly spaced values within a specified range. In this case, from 0 to TDur in steps of dt
Nt = len(t)  # the length of t

Eta = np.zeros(Nt)
U = np.zeros((Nz, Nt))  # horizontal velocity
Ut = np.zeros((Nz, Nt))  # horizontal acceleration

Fx = np.zeros(Nt)  # force in x
My = np.zeros(Nt)  # moment in y

for it in range(0, Nt):
    # in the irregular wave formula it should be like "a*np.cos(omega*t[it]-k*x+epsilon)", but we choose x to be 0 because we are only looking at the surface of the monopile
    # this is equivalent as using the regular wave formula for this specific conditions
    Eta[it] = sum(a * np.cos(omega * t[it] + waves['epsilon']));
    # zPhys is an evenly spaced array that goes from the bottom of the seabed (-h) to the surface elevation (Eta[it]). The size of the array is Nz
    zPhys = np.linspace(-h, Eta[it], Nz);
    zCalc = np.linspace(-h, 0., Nz) # Uuse zCalc in the calculation of velocity and acceleration to apply Wheeler stretching

    for iz in range(Nz):
        # Calculates the horizontal velocity
        U[iz, it] = np.sum(a * omega * np.cosh(k * (zCalc[iz] + h)) / np.sinh(k * h) * np.cos(omega * t[it] + waves['epsilon']))
        Ut[iz, it] = np.sum(-1*a * omega**2 * np.cosh(k * (zCalc[iz] + h)) / np.sinh(k * h) * np.sin(omega * t[it] + waves['epsilon']))

    [F, M] = forceIntegrate(monopileData, U[:, it], Ut[:, it], zPhys)

    Fx[it] = F
    My[it] = M

# Make some statistics
print("{:8} {:^15} {:^15} {:^15}".format('', 'eta', 'Fx', 'My'))
columnFormat = "{:8} {:^15.2f} {:^15.2f} {:^15.2f}"
print(columnFormat.format('mean', np.mean(Eta), np.mean(Fx), np.mean(My)))
print(columnFormat.format('std', np.std(Eta), np.std(Fx), np.std(My)))
print(columnFormat.format('min', np.min(Eta), np.min(Fx), np.min(My)))
print(columnFormat.format('max', np.max(Eta), np.max(Fx), np.max(My)))

print("4 sigma_eta:", 4*np.std(Eta))
print("Hs",Hs)
f2, ax2 = plt.subplots(3, sharex=True)

ax2[0].plot(t, Eta, 'b')
ax2[0].set_xlabel('t [s]')
ax2[0].set_ylabel('η[m]')
ax2[0].grid()

ax2[1].plot(t, Fx/10e6, 'k')
ax2[1].set_xlabel('t [s]')
ax2[1].set_ylabel('Force [N]')
ax2[1].grid()

ax2[2].plot(t, My, 'r')
ax2[2].set_xlabel('t [s]')
ax2[2].set_ylabel('Moment [Nm]')
ax2[2].grid()

f2.tight_layout()


plt.figure()
plt.plot(t, Fx/10e6)
plt.xlabel("Time [s]")
plt.ylabel("Force [MN]")
plt.xlim(0, TDur)
plt.grid(True)
plt.show()

np.savetxt('Myassignment1.txt', My)
np.savetxt('MYASSJS.txt', S)

#x1, y1 = get_plot1_data()
#np.savez('plot1_data.npz', x1=x1, y1=y1)