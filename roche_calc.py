import numpy as np

# Global Constants
G = 6.67e-11 # G in SI units ****
M_SUN = 1.9891e30 # mass of sun in kg ****
AU_TO_M = 1.496e+11 # meters in one AU
YR_TO_SEC = 365.25*24.*3600. # seconds in a year

rSat = 60268000. # meters

rhoSat = 1.88*1000. # kg/m^3
rhoTitan = 687 # kg/m^3

roche_r = 1.26*rSat*np.power(rhoSat/rhoTitan, 1./3.)

a = 12.*rSat

print(1-(roche_r/a))