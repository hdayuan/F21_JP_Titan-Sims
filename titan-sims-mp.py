"""
titan-sims-mp.py
Author: Henry Yuan

Description:
Integrates the Iapetus-Titan-Saturn-Sun system as Titan migrates
through the evection resonance. Includes Saturn's obliquity and the tidal
effects on Saturn and Titan, but not on the Sun. Includes Iapetus, to see
if its inclination is affected.

To run:

python3 titan-sims-mp.py [continuation?] [initial a] [number of samples]
[total integration time] {[number of additional samples]
[additional integration time]}

Takes 4-6 command line arguments: continuation is 0 if new sim,
1 if continuing previous sim, number of samples, initial semi-major axis of 
Titan in Saturn radii, and total simulation time in millions of years, {number
of additional samples, additional integration time in millions of years}

*** Units are AU, yr, and MSun (and SI for other measurements) 
unless otherwise specified ***

Trials
"""

import multiprocessing as mp
import sys
import os
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


"""Creates a new file for a new simulation, writes the parameters of simulation
to the file, and returns the file"""
def new_sim_file(iaTitanRS, numSamples, intTime, k2Titan):
    file_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs-k"+str(k2Titan)
    
    f = open("v4.4-"+file_str+".txt", "a")

    # Write parameters of simulation
    f.write(str(iaTitanRS)+" Saturn radii\n")
    f.write(str(numSamples)+" samples\n")
    f.write(str(intTime)+" million years\n")
    f.write("Titan's k2 = "+str(k2Titan)+"\n")

    return f


"""Integrates the system starting from a = iaTitanRS saturn radii for intTime 
million years and prints numSamples data points to file, 
each on one line in the following format:

'a (RS)'[tab]'e'[tab]'i (rad)'[tab]'longitude of pericenter'[tab]
'mean longitude of Sun'[tab]'Iapetus a'
[tab]'Iapetus e'[tab]'Iapetus i'[tab]'current time in simulation'

Lastly, saves simulation"""
def integrate_sim(iaTitanRS, numSamples, intTime, k2Titan, file):

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
    eTitan = 0.001 # initial eccentricity
    rTitan = 0.04421567543 * rSat # radius
    QTitan = 100. # Estimated tidal Q factor (Murray & Dermott p 173)

    # Iapetus constants
    mIap = 1.806e21 / M_SUN # mass
    aIap = 3561000000. / AU_TO_M # modern-day semi-major axis

    # more constants
    iaTitan = iaTitanRS * rSat  # starting semi-major axis of Titan
    aResRS = evection_a(mSat, rSat, j2Sat, aSat) # Titan's semi-major axis at resonance in units of saturn radii
    nTitan = mean_motion(mSat, aResRS*rSat) # mean motion of Titan at resonance in rad / sec
    omegaTitan = nTitan # spin rate of Titan, in rad/sec (same as nTitan because Titan is tidally locked)
    perTitan = get_period(nTitan) # period of Titan at resonance
    timescale = get_timescale(aResRS*rSat, aTitan, ageSat) # migration timescale at resonance

    t_start = ageSat * (iaTitan/aTitan)**3
    iaIap = aIap * (t_start/ageSat)**(1./3.) # initial semi-major axis of Iapetus

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    sim.dt = (1./20.) * perTitan # time step = 1/20 * shortest orbital period

    # add Saturn, Titan, Iapetus, and Sun
    saturn = rebound.Particle(m=mSat, hash='Saturn')
    sim.add(saturn)
    titan = rebound.Particle(sim, primary=saturn, m=mTitan, a=iaTitan, e=eTitan, hash='Titan')
    sim.add(titan)
    iapetus = rebound.Particle(sim, primary=saturn, m=mIap, a=iaIap, e=0, hash='Iapetus')
    sim.add(iapetus)
    sun = rebound.Particle(sim, primary=saturn, m=1., a=aSat, e=eSat, inc=oSat, hash='Sun')
    sim.add(sun)

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
    ps['Titan'].params["tctl_tau"] = time_lag_tau(nTitan, 0., QTitan) # yrs
    ps['Titan'].params["Omega"] = omegaTitan*YR_TO_SEC # rad/yr

    # Tidal forces of Saturn
    QSat = 3.*k2Sat*(mTitan/mSat)*nTitan*(timescale*YR_TO_SEC)*(1./aResRS)**5. # Q for Saturn based on timescale
    ps['Saturn'].r = rSat # AU
    ps['Saturn'].params["tctl_k2"] = k2Sat
    ps['Saturn'].params["tctl_tau"] = time_lag_tau(nTitan, omegaSat, QSat) # yrs
    ps['Saturn'].params["Omega"] = omegaSat*YR_TO_SEC # rad/yr

    # Integrate
    plotDT = intTime*1000000./numSamples
    for i in range(numSamples):
        sim.integrate(i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(ps['Titan'].a/rSat)+"\t"+str(ps['Titan'].e)+"\t"+str(ps['Titan'].inc)+"\t")
        file.write(str(ps['Titan'].pomega)+"\t"+str(ps['Sun'].l)+"\t")
        file.write(str(ps['Iapetus'].a/rSat)+"\t"+str(ps['Iapetus'].e)+"\t")
        file.write(str(ps['Iapetus'].inc)+"\t")
        file.write(str(sim.t)+"\n")

    # save simulation
    # file_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
    # sim.save("v4.4-sim-"+file_str+".bin")
    # rebx.save("v4.4-simx-"+file_str+".bin")

    # sim.cite()


"""Opens file, writes command line arguments, starts timer, 
calls integration method, stops timer and writes real-time duration
of simulation"""
def run_sim(iaTitanRS, numSamples, intTime, trial):

    k2Titan = trial * 0.05
    f = new_sim_file(iaTitanRS, numSamples, intTime, k2Titan)

    # start timer
    start_time = time.time()



    integrate_sim(iaTitanRS, numSamples, intTime, k2Titan, f)
    
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
# Step 1: Init multiprocessing.Pool()
numCPUs = mp.cpu_count()
#pool = mp.Pool(numCPUs)

# Step 2: `pool.apply` the `run_sim()`
#results = [pool.apply(run_sim, args=(8.2, 4000, 20.0, trial)) for trial in numCPUs]

# Step 3: Don't forget to close
#pool.close()    

#print(results[:10])