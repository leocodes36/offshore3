import numpy as np

def forceIntegrate(monopileData, u, ut, zPhys):
    
    df = 0.5*monopileData["rho"]*monopileData["D"]*monopileData["CD"]*np.abs(u)*u + monopileData["rho"]*monopileData["D"]*monopileData["D"]*np.pi/4*monopileData["CM"]*ut
    F = np.trapz(df, zPhys)
    M = np.sum(df*(monopileData["h"]+zPhys)) # moment is defined as sum of all moments (each force with its respective lever), h of monopile is added to zphys to make lever positive becuase moment acts in the same direction as force (zphys has negative values due to position of z at swl)

    return F, M


