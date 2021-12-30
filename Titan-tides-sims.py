"""Description:
Integrates the Titan-Saturn-Sun system as Titan migrates
through the evection resonance. Includes Saturn's obliquity and the tidal
effects on Saturn and Titan, but not on the Sun. Also includes Jupiter which
causes oscillations in Saturn's eccentricity

TO-DO:
-print exp a res
-Check values of k2 for Saturn and Titan
-how to set inclination of Jupiter?
-use different equatorial radius for Saturn?
-figure out how to save simulations so that I can pick up where
the last one left off
-figure out how sim.cite() works and whether to use it
-see JPNotes doc and JP_introduction doc for other questions/to-dos

To run:

python3 Titan-sims.py [number of samples] [initial a] [total integration time]

Takes 3 command line arguments: number of samples, initial semi-major axis of 
Titan in Saturn radii, and total simulation time in in millions of years

*** Units are AU, yr, and MSun (and SI for other measurements) 
unless otherwise specified ***"""


import sys
import rebound
import reboundx
import numpy as np
import time


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


"""Integrates the system starting from a = iaTitanRS saturn radii for intTime 
years and prints numSamples data points to file, 
each on one line in the following format:

'a (RS)'[tab]'e'[tab]'i (rad)'[tab]'longitude of pericenter'[tab]
'mean longitude of Sun'[tab]'e of Sun'[tab]'current time in simulation'

Lastly, saves simulation"""
def integrate_sim(numSamples, iaTitanRS, intTime, file):

    # Check accuracy of all the following constants
    # Saturn constants
    mSat = 0.0002857 # mass
    aSat = 9.5549 # semi-major axis
    eSat = 0.0565 # eccentricity
    # reqSat = 60268000./AU_TO_M # equatorial radius--need this?
    rSat = 58232503./AU_TO_M # mean physical radius
    oSat = 26.7 * np.pi / 180. # obliquity
    j2Sat = 16298e-6 # J2 (Murray and Dermott p 531)
    k2Sat = 0.341 # Love number *** CHECK THIS ***
    omegaSat = 2*np.pi/(10.656*3600) # spin rate (*** rad/sec ***)
    ageSat = 4.503e9 # age

    # Jupiter Constants
    # mJup = (1.898e+27)/M_SUN
    # aJup = (778.570e+9)/AU_TO_M
    # eJup = 0.0489

    # Titan constants
    mTitan = 0.0000000676319759 # mass
    aTitan = 0.008167696467 # modern-day semi-major axis
    eTitan = 0.001 # initial eccentricity
    rTitan = 0.04421567543 * rSat # radius
    k2Titan = 0.15 # Love number *** CHECK THIS ***
    QTitan = 100. # Estimated tidal Q factor

    # more constants
    iaTitan = iaTitanRS * rSat  # starting semi-major axis of Titan
    aResRS = evection_a(mSat, rSat, j2Sat, aSat) # Titan's semi-major axis at resonance in units of saturn radii
    nTitan = mean_motion(mSat, aResRS*rSat) # mean motion of Titan at resonance in rad / sec
    omegaTitan = nTitan # spin rate of Titan, in rad/sec (same as nTitan because Titan is tidally locked)
    perTitan = get_period(nTitan) # period of Titan at resonance
    timescale = get_timescale(aResRS*rSat, aTitan, ageSat) # migration timescale at resonance

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    sim.dt = (1./20.) * perTitan # initial time step = 1/20 * shortest orbital period

    # add Saturn, Titan, Sun, and Jupiter (what inclination?)
    saturn = rebound.Particle(m=mSat, hash='Saturn')
    sim.add(saturn)
    titan = rebound.Particle(sim, primary=saturn, m=mTitan, a=iaTitan, e=eTitan, hash='Titan')
    sim.add(titan)
    sun = rebound.Particle(sim, primary=saturn, m=1., a=aSat, e=eSat, inc=oSat, hash='Sun')
    sim.add(sun)
    # jupiter = rebound.Particle(sim, primary=sun, m=mJup, a=aJup, e=eJup, hash='Jupiter')
    # sim.add(jupiter)

    # Initiate reboundx
    rebx = reboundx.Extras(sim)
    ps = sim.particles

    # add Saturn's J2
    gh = rebx.load_force("gravitational_harmonics")
    rebx.add_force(gh)
    ps['Saturn'].params["J2"] = j2Sat
    ps['Saturn'].params["R_eq"] = rSat

    # add tidal forces
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    # Tidal forces of Titan
    ps['Titan'].r = rTitan # AU
    ps['Titan'].params["tctl_k2"] = k2Titan
    ps['Titan'].params["tctl_tau"] = time_lag_tau(nTitan, omegaTitan, QTitan) # yrs
    ps['Titan'].params["Omega"] = omegaTitan*YR_TO_SEC # rad/yr

    # Tidal forces of Saturn
    QSat = 3.*k2Sat*(mTitan/mSat)*nTitan*(timescale*YR_TO_SEC)*(1./aResRS)**5. # Q for Saturn based on timescale
    ps['Saturn'].r = rSat # AU
    ps['Saturn'].params["tctl_k2"] = k2Sat
    ps['Saturn'].params["tctl_tau"] = time_lag_tau(nTitan, omegaSat, QSat) # yrs
    ps['Saturn'].params["Omega"] = omegaSat*YR_TO_SEC # rad/yr

    # Integrate
    plotDT = intTime/numSamples
    for i in range(numSamples):
        sim.integrate(i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(ps['Titan'].a/rSat)+"\t"+str(ps['Titan'].e)+"\t"+str(ps['Titan'].inc)+"\t")
        # file.write(str(sun.e)+"\t")
        file.write(str(ps['Titan'].pomega)+"\t"+str(ps['Sun'].l)+"\t")
        file.write(str(sim.t)+"\n")

    # sim.cite()
    # save simulation


"""Opens file, writes command line arguments, starts timer, 
calls integration method, stops timer and writes real-time duration
of simulation"""
def main():

    # Parse command-line arguments
    numSamples = int(sys.argv[1])
    iaTitanRS = float(sys.argv[2])
    intTime = float(sys.argv[3]) # total time to integrate, in myrs

    # open file
    f = open("v4t-"+str(numSamples)+"s-"+str(iaTitanRS)+"rs-"+str(intTime)+"myrs.txt", "a")

    # Write parameters of simulation
    f.write(str(numSamples)+" samples\n")
    f.write(str(iaTitanRS)+" Saturn radii\n")
    f.write(str(intTime)+" million years\n")

    # start timer
    start_time = time.time()

    # call integration method
    integrate_sim(numSamples, iaTitanRS, intTime*1000000., f)

    # Write running time
    totTimeSec = time.time() - start_time
    numDays = int(totTimeSec) // (3600*24)
    numHours = int(totTimeSec) // 3600
    numMins = int(totTimeSec % 3600) // 60
    numSecs = (totTimeSec % 3600) % 60
    f.write("Running time: "+str(numDays)+ " days "+str(numHours)+ " hours "+str(numMins)+" minutes "+str(numSecs)+" seconds.\n")

    # close the file
    f.close()

# Call main()
main()
