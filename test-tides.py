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


"""Calculates the mean motion (rad/sec) of a body orbiting around a central body
of mass m (MSun) at semi-major axis a (AU)"""
def mean_motion(m, a):
    return np.sqrt(G*m*M_SUN / ((a*AU_TO_M)**3))


"""Calculates the period (yrs) of a satellite orbiting a body
with mean motion n (rad/sec)"""
def get_period(n):
    return 2*np.pi/n/YR_TO_SEC


"""Calculates the semi-major axis at which evection resonance occurs for 
a satellite orbiting a planet of mass mp (MSun), radius rp (AU), and J2 moment j2p,
that is itself orbiting the Sun at a semi-major axis a (AU).
Assumes i = e = 0. Return value in units of the planet's radius"""
def evection_a(mp, rp, j2p, a):
    return np.power((9./4.)*j2p**2*mp*(a/rp)**3, (1./7.))


"""Calculates timescale (yr) for migration of semi-major axis for satellite at 
semi-major axis a (AU) but whose modern-day semi-major axis is aCurr (AU),
orbiting planet whose age is age (yr)"""
def get_timescale(a, aCurr, age):
    return 3. * ((a/aCurr)**3) * age


"""COMMENT"""
def main(ia_titanRS, simTime, eTitan, incTitan):

    # Check accuracy of all the following constants
    # Saturn constants
    mSat = 0.0002857 # mass
    aSat = 9.5549 # semi-major axis
    eSat = 0.0565 # eccentricity
    rSat = 58232503./AU_TO_M # mean physical radius
    oSat = 26.7 * np.pi / 180. # obliquity
    j2Sat = 16298e-6 # J2 (Murray and Dermott p 531)
    k2Sat = 0.341 # from Cuk et al. 2016
    omegaSat = 2*np.pi/(10.656*3600) # spin rate (*** rad/sec ***)
    ageSat = 4.503e9 # age

    # Titan constants
    mTitan = 0.0000000676319759 # mass
    aTitan = 0.008167696467 # modern-day semi-major axis
    rTitan = 0.04421567543 * rSat # radius
    k2Titan = 0.15 # Love number (calculated from Murray & Dermott p 173)
    QTitan = 100. # Estimated tidal Q factor (Murray & Dermott p 173)

    # more constants
    iaTitan = ia_titanRS * rSat  # starting semi-major axis of Titan
    aResRS = evection_a(mSat, rSat, j2Sat, aSat) # Titan's semi-major axis at resonance in units of saturn radii
    nTitan = mean_motion(mSat, aResRS*rSat) # mean motion of Titan at resonance in rad / sec
    omegaTitan = nTitan # spin rate of Titan, in rad/sec (same as nTitan because Titan is tidally locked)
    perTitan = get_period(nTitan) # period of Titan at resonance
    timescale = get_timescale(aResRS*rSat, aTitan, ageSat) # migration timescale at resonance

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    nTitan = mean_motion(mSat, iaTitan)
    tauTitan = 2.*np.pi/(nTitan*YR_TO_SEC)
    sim.dt = (1./20.) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat, hash = "Saturn")

    # add Titan
    sim.add(m=mTitan, a=iaTitan, e=eTitan, inc=incTitan, hash="Titan")

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

    # Tidal forces of Saturn
    QSat = 3.*k2Sat*(mTitan/mSat)*nTitan*(timescale*YR_TO_SEC)*(1./aResRS)**5. # Q for Saturn based on timescale
    ps['Saturn'].r = rSat # AU
    ps['Saturn'].params["tctl_k2"] = k2Sat
    ps['Saturn'].params["tctl_tau"] = time_lag_tau(nTitan, omegaSat, QSat) # yrs
    ps['Saturn'].params["Omega"] = omegaSat*YR_TO_SEC # rad/yr

    sim.integrate(simTime*1000000.)
    # move to Saturn's frame of reference
    sim.move_to_hel()
    # write output
    print(ps["Titan"].inc*180/np.pi)

    # inclination damping calculations
    mass_ratio = 4225.
    mS = mSat*M_SUN
    rS = rSat*AU_TO_M
    Omega = omegaSat

    i0 = incTitan
    a0 = 8.2*rS

    const_1 = np.sqrt(1+(1/mass_ratio))/4.
    const_2 = np.sqrt(1+(1/mass_ratio))/6.*np.sqrt(G*mS)/Omega

    # calculate last constant based on data point above
    const = i0*a0**const_1/(np.exp(const_2*a0**(-3./2.)))

    a = ps['Titan'].a*AU_TO_M
    i = const/(a**const_1)*np.exp(const_2*a**(-3./2.)) # radians
    i_deg = i*180./np.pi
    print(i_deg)

    print(ps['Titan'].e)

    """numSteps = simTime*1000000./sim.dt
    e = eTitan
    e_lim = e
    a = iaTitan*AU_TO_M
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
    print(e_lim)"""

# Parse command-line arguments
simTime = float(sys.argv[1]) # myr
eTitan = float(sys.argv[2])
incTitan = float(sys.argv[3]) # degrees

main(8.2, simTime, eTitan, incTitan*np.pi/180.)
