import sys
import rebound
import reboundx
import numpy as np

# Global Constants
G = 6.67e-11 # G in SI units ****
M_SUN = 1.9891e30 # mass of sun in kg ****
AU_TO_M = 1.496e+11 # meters in one AU
YR_TO_SEC = 365.25*24.*3600. # seconds in a year

"""Calculates tidal time-lag T (in years) for a body with mean motion n (radians/sec), 
tidal Q factor Q, radius r (AU), and mass m (solar masses)"""
def time_lag(n, Q, r, m):
    return 2*n*Q*(r*AU_TO_M)**3 / (G*m*M_SUN*YR_TO_SEC)

"""Calculates the mean motion (in radians/sec) of a body orbiting around a central body
of mass mCentral (solar masses) at semi-major axis a (AU) """
def mean_motion(mCentral, a):
    return np.sqrt(G*mCentral*M_SUN / ((a*AU_TO_M)**3))

"""COMMENT"""
def main(ia_titanRS, fa_titanRS):

    # Check accuracy of all the following constants
    mSat = 0.0002857 # mass of saturn in solar masses
    rSat = 0.00038926024 # radius of Saturn in AU
    omegaSat = 2*np.pi*YR_TO_SEC/(10.656*3600) # spin rate of Saturn in radians per year
    j2Sat = 16298e-6 # J2 of Saturn (Murray and Dermott p 531)
    # calculated later # QSat = 5000. # Tidal Q factor of Saturn (Lainey et al.)
    k2Sat = 0.39 # Love number of Saturn *** CHECK THIS ***
    curr_aTitan = 0.008167696467 # modern-day semi-major axis of titan in AU
    mTitan = 0.0000000676319759 # mass of titan in solar masses
    rTitan = 0.04421567543 * rSat # radius of Titan in AU
    k2Titan = 0.15 # Love number of Titan *** CHECK THIS ***
    QTitan = 100. # Tidal Q factor of Titan

    ia_titanAU = ia_titanRS * rSat  # starting semi-major axis of Titan

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    nTitan = np.sqrt(G*mSat*M_SUN / ((8.215*rSat*AU_TO_M)**3))
    tauTitan = 2.*np.pi/(nTitan*YR_TO_SEC)
    sim.dt = (1./20.) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat, hash = "Saturn")

    # add Titan
    sim.add(m=mTitan, a=ia_titanAU, e=0, hash = "Titan")

    timescale = 811000000.

    # calculate total time to integrate based on ia_titan, fa_titan, tau, and
    # given exponential migration of form a = a_0 * e^(t/tau)
    totSimTime = timescale * np.log(fa_titanRS / ia_titanRS)

    # Initiate reboundx
    rebx = reboundx.Extras(sim)
    ps = sim.particles

    # add Saturn's J2
    # gh = rebx.load_force("gravitational_harmonics")
    # rebx.add_force(gh)
    # ps["Saturn"].params["J2"] = j2Sat
    # ps["Saturn"].params["R_eq"] = rSat

    # add tidal forces
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    # Tidal forces of Titan
    # titanT = time_lag(nTitan, QTitan, rTitan, mTitan) # in years
    # ps["Titan"].r = rTitan # AU
    # ps["Titan"].params["tctl_k2"] = k2Titan
    # ps["Titan"].params["tctl_tau"] = titanT
    # ps["Titan"].params["Omega"] = nTitan*YR_TO_SEC # in radians per year

    # Tidal forces of Saturn
    QSat = 3.*k2Sat*(mTitan/mSat)*(omegaSat/YR_TO_SEC)*timescale*(1.0/8.215)**5. # Q for Saturn 
    satT = time_lag((omegaSat/YR_TO_SEC), QSat, rSat, mSat) # in years
    ps["Saturn"].r = rSat # AU
    ps["Saturn"].params["tctl_k2"] = k2Sat
    ps["Saturn"].params["tctl_tau"] = satT
    ps["Saturn"].params["Omega"] = omegaSat # in rad per year

    sim.integrate(totSimTime)
    # move to Saturn's frame of reference
    sim.move_to_hel()
    # write output
    print(str(ps["Titan"].a/rSat))

# Parse command-line arguments
ia_titan = float(sys.argv[1])
fa_titan = float(sys.argv[2])

main(ia_titan, fa_titan)
