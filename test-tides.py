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
    return 1./(2.*Q*np.abs(omega-n))/YR_TO_SEC

"""Calculates the mean motion (in radians/sec) of a body orbiting around a central body
of mass mCentral (solar masses) at semi-major axis a (AU) """
def mean_motion(mCentral, a):
    return np.sqrt(G*mCentral*M_SUN / ((a*AU_TO_M)**3))

"""COMMENT"""
def main(ia_titanRS, simTime, eTitan):

    # Check accuracy of all the following constants
    mSat = 0.0002857 # mass of saturn in solar masses
    rSat = 0.00038926024 # radius of Saturn in AU
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
    sim.add(m=mTitan, a=ia_titanAU, e=eTitan, hash="Titan")

    # Initiate reboundx
    rebx = reboundx.Extras(sim)
    ps = sim.particles

    # add tidal forces
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    # Tidal forces of Titan
    ps['Titan'].r = rTitan # AU
    ps['Titan'].params["tctl_k2"] = k2Titan
    ps['Titan'].params["tctl_tau"] = time_lag_tau(nTitan, 0., QTitan) # yrs
    ps['Titan'].params["Omega"] = nTitan*YR_TO_SEC # rad/yr

    sim.integrate(simTime*1000000.)
    # move to Saturn's frame of reference
    sim.move_to_hel()
    # write output
    print(str(ps["Titan"].e))

    numSteps = simTime*1000000./sim.dt
    e = eTitan
    e_lim = e
    a = ia_titanAU*AU_TO_M
    for i in range(int(numSteps)):
        # Calculate expected e damp (assuming linear damping, for small
        # time intervals)
        q = mSat/mTitan
        R_s = rTitan*AU_TO_M
        M_s = mTitan*M_SUN
        M_p = mSat*M_SUN
        T = R_s**3/(G*M_s)*2.*nTitan*QTitan
        f3 = 1 + (15./4.)*(e**2) + (15./8.)*(e**4) + (5./64.)*(e**6)
        f4 = 1 + (3./2.)*(e**2) + (1./8.)*(e**4)
        e_dot = -27.*k2Titan/T*q*(1.+q)*(R_s/a)**8*e/((1.-e**2)**(13./2.))*(f3-(11./18.)*(1.-e**2)**(3./2.)*f4)
        e = e+(e_dot*sim.dt*YR_TO_SEC)

        # calculate small e limit expected final e
        e_dot_lim = -(63./8.)*e_lim*nTitan/900.*M_p/M_s*(R_s/a)**5
        e_lim = e_lim+(e_dot_lim*sim.dt*YR_TO_SEC)

    print(e)
    print(e_lim)
# Parse command-line arguments
simTime = float(sys.argv[1]) # myr
eTitan = float(sys.argv[2])

main(8.2, simTime, eTitan)
