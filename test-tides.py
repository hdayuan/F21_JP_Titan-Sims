import sys
import rebound
import reboundx
import numpy as np

# Global Constants
G = 6.67e-11 # G in SI units ****
M_SUN = 1.9891e30 # mass of sun in kg ****
AU_TO_M = 1.496e+11 # meters in one AU
YR_TO_SEC = 365.25*24.*3600. # seconds in a year

"""Calculates tau (yrs) for constant time lag for a body with spin rate omega
(rad/sec), a satellite orbiting with mean motion n (rad/sec), tidal Q factor Q"""
def time_lag_tau(n, omega, Q):
    if n == omega:
        return 0
    return 1./(2.*Q*np.abs(omega-n))/YR_TO_SEC

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
    k2Sat = 0.341 # Love number of Saturn *** CHECK THIS ***
    mTitan = 0.0000000676319759 # mass of titan in solar masses
    rTitan = 0.04421567543 * rSat # radius of Titan in AU
    k2Titan = 0.15 # Love number of Titan *** CHECK THIS ***
    QTitan = 100. # Tidal Q factor of Titan

    ia_titanAU = ia_titanRS * rSat  # starting semi-major axis of Titan

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    nTitan = mean_motion(mSat, ia_titanAU)
    tauTitan = 2.*np.pi/(nTitan*YR_TO_SEC)
    sim.dt = (1./20.) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat, hash = "Saturn")

    # add Titan
    sim.add(m=mTitan, a=ia_titanAU, e=0, hash = "Titan")

    timescale = 810766627.6731131 # in years

    # calculate total time to integrate based on ia_titan, fa_titan, tau, and
    # given exponential migration of form a = a_0 * e^(t/tau)
    totSimTime = timescale * np.log(fa_titanRS / ia_titanRS)

    # Initiate reboundx
    rebx = reboundx.Extras(sim)
    ps = sim.particles

    # add tidal forces
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    # Tidal forces of Saturn
    QSat = 3.*k2Sat*(mTitan/mSat)*nTitan*(timescale*YR_TO_SEC)*(1.0/8.215)**5. # Q for Saturn 

    ps["Saturn"].r = rSat # AU
    ps["Saturn"].params["tctl_k2"] = k2Sat
    ps["Saturn"].params["tctl_tau"] = time_lag_tau(nTitan, omegaSat/YR_TO_SEC, QSat)
    ps["Saturn"].params["Omega"] = omegaSat #/4.4013 # in rad per year

    sim.integrate(totSimTime)
    # move to Saturn's frame of reference
    sim.move_to_hel()
    # write output
    print(str(ps["Titan"].a/rSat))

    # Calculate expected migration distance (assuming linear migration, for small
    # time intervals)

    # calculate migration rate at initial a
    a_dot = 3.*k2Sat*mTitan/(QSat*mSat)*(1./ia_titanRS)**5*nTitan*ia_titanAU*AU_TO_M
    final_a = (ia_titanAU+(a_dot*totSimTime*YR_TO_SEC)/AU_TO_M)/rSat
    print(final_a)

# Parse command-line arguments
ia_titan = float(sys.argv[1])
fa_titan = float(sys.argv[2])

main(ia_titan, fa_titan)
