"""
titan-sims-v4.2.py
Author: Henry Yuan

Description:
Integrates the Titan-Saturn-Sun system as Titan migrates
through the evection resonance. Includes Saturn's obliquity and the tidal
effects on Saturn and Titan, but not on the Sun.

To run:

python3 Titan-tides-sims.py [continuation?] [initial a] [number of samples]
[total integration time] {[number of additional samples]
[additional integration time]}

Takes 4-6 command line arguments: continuation is 0 if new sim,
1 if continuing previous sim, number of samples, initial semi-major axis of 
Titan in Saturn radii, and total simulation time in millions of years, {number
of additional samples, additional integration time in millions of years}

*** Units are AU, yr, and MSun (and SI for other measurements) 
unless otherwise specified ***

TO-DO:
-Check values of k2 for Saturn and Titan
-how to set inclination of Jupiter?
-use different equatorial radius for Saturn?
-figure out how sim.cite() works and whether to use it
-see JPNotes doc and JP_introduction doc for other questions/to-dos
"""


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


"""Creates a new file for a new simulation, writes the parameters of simulation
to the file, and returns the file"""
def new_sim_file(iaTitanRS, numSamples, intTime):
    file_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
    
    f = open("v4.2-"+file_str+".txt", "a")

    # Write parameters of simulation
    f.write(str(iaTitanRS)+" Saturn radii\n")
    f.write(str(numSamples)+" samples\n")
    f.write(str(intTime)+" million years\n")
    f.write("Continuation of: none\n")

    return f


"""Creates a new file for a continued simulation, writes the parameters of simulation
to the file, copies data from previous sim file, and returns the file"""
def cont_sim_file(iaTitanRS, numSamples, intTime, numSamplesADD, intTimeADD):
    file_prev_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
    file_new_str = str(iaTitanRS)+"rs-"+str(numSamples+numSamplesADD)+"s-"+str(intTime+intTimeADD)+"myrs"

    f_prev = open(r"v4.2-"+file_prev_str+".txt", "r")
    f_new = open("v4.2-"+file_new_str+".txt", "a")
    
    # Write parameters of simulation
    f_new.write(str(iaTitanRS)+" Saturn radii\n")
    f_new.write(str(numSamples+numSamplesADD)+" samples\n")
    f_new.write(str(intTime+intTimeADD)+" million years\n")
    f_new.write("Continuation of: v4.2-"+file_prev_str+".txt\n")

    # copy data from f_prev
    f_prev.readline()
    f_prev.readline()
    f_prev.readline()
    f_prev.readline()

    for i in range(numSamples):
        f_new.write(f_prev.readline())

    return f_new


"""Continues a previous simulation whose parameters are iaTitanRS, numSamples,
and intTime, integrating for intTimeADD million more years and taking numSamplesADD
more samples, appending new data to file. See integrate_sim description for
format of output. Lastly, saves simulation"""
def continue_sim(iaTitanRS, numSamples, intTime, numSamplesADD, intTimeADD, file):
    
    rSat = 58232503./AU_TO_M # mean physical radius
    file_prev_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
    file_new_str = str(iaTitanRS)+"rs-"+str(numSamples+numSamplesADD)+"s-"+str(intTime+intTimeADD)+"myrs"
    
    # reload previous simulation
    sim = rebound.Simulation("v4.1-sim-"+file_prev_str+".bin")
    rebx = reboundx.Extras(sim, "v4.1-simx-"+file_prev_str+".bin")
    ps = sim.particles
    sim_start_time = sim.t # current time in previous simulation

    # Integrate
    plotDT = intTimeADD*1000000./numSamplesADD
    for i in range(1, numSamplesADD + 1):
        sim.integrate(sim_start_time + i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(ps['Titan'].a/rSat)+"\t"+str(ps['Titan'].e)+"\t"+str(ps['Titan'].inc)+"\t")
        # file.write(str(sun.e)+"\t")
        file.write(str(ps['Titan'].pomega)+"\t"+str(ps['Sun'].l)+"\t")
        file.write(str(sim.t)+"\n")

    # remove old binary files of saved simulation
    os.remove("v4.1-sim-"+file_prev_str+".bin")
    os.remove("v4.1-simx-"+file_prev_str+".bin")

    # save simulation
    sim.save("v4.1-sim-"+file_new_str+".bin")
    rebx.save("v4.1-simx-"+file_new_str+".bin")

    # sim.cite()


"""Integrates the system starting from a = iaTitanRS saturn radii for intTime 
million years and prints numSamples data points to file, 
each on one line in the following format:

'a (RS)'[tab]'e'[tab]'i (rad)'[tab]'longitude of pericenter'[tab]
'mean longitude of Sun'[tab]'e of Sun'[tab]'current time in simulation'

Lastly, saves simulation"""
def integrate_sim(iaTitanRS, numSamples, intTime, file):

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
    plotDT = intTime*1000000./numSamples
    for i in range(numSamples):
        sim.integrate(i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(ps['Titan'].a/rSat)+"\t"+str(ps['Titan'].e)+"\t"+str(ps['Titan'].inc)+"\t")
        # file.write(str(sun.e)+"\t")
        file.write(str(ps['Titan'].pomega)+"\t"+str(ps['Sun'].l)+"\t")
        file.write(str(sim.t)+"\n")

    # save simulation
    file_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
    sim.save("v4.1-sim-"+file_str+".bin")
    rebx.save("v4.1-simx-"+file_str+".bin")

    # sim.cite()


"""Opens file, writes command line arguments, starts timer, 
calls integration method, stops timer and writes real-time duration
of simulation"""
def main():

    # Parse command-line arguments
    continuation = int(sys.argv[1]) # 0 if new sim, 1 if continuing a previous sim
    iaTitanRS = float(sys.argv[2])
    numSamples = int(sys.argv[3])
    intTime = float(sys.argv[4]) # total time to integrate, in myrs

    # Check validity of command-line arguments
    if (continuation != 0 and continuation != 1):
        raise Exception("First command-line argument must be 0 or 1: 0 if it is "
        +"a new simulation, 1 if it is a continuation of a previous simuation")
    if (continuation == 1):
        file_str = str(iaTitanRS)+"rs-"+str(numSamples)+"s-"+str(intTime)+"myrs"
        if (not os.path.exists("v4.1-sim-"+file_str+".bin")) or (not os.path.exists("v4.1-simx-"+file_str+".bin")):
            raise Exception("One or more binary files for previous simulation to be continued does not exist")
        if not os.path.exists("v4.2-"+file_str+".txt"):
            raise Exception("Data file for previous simulation to be continued does not exist")

    # if continuation, 2 more command-line args
    if (continuation == 1):
        numSamplesADD = int(sys.argv[5]) # additional number of samples
        intTimeADD = float(sys.argv[6]) # additional time to integrate, in myrs
        # if numSamplesADD is 0, use same ratio of samples as previous simulation
        if (numSamplesADD == 0):
            numSamplesADD = int(intTimeADD*numSamples/intTime)

    # open file
    if (continuation == 0):
        f = new_sim_file(iaTitanRS, numSamples, intTime)
    else:
        f = cont_sim_file(iaTitanRS, numSamples, intTime, numSamplesADD, intTimeADD)

    # start timer
    start_time = time.time()

    # call integration method
    if (continuation == 0):
        integrate_sim(iaTitanRS, numSamples, intTime, f)
    else:
        continue_sim(iaTitanRS, numSamples, intTime, numSamplesADD, intTimeADD, f)

    # Write running time (if continuation, only running time for added integration)
    totTimeSec = time.time() - start_time
    numDays = int(totTimeSec) // (3600*24)
    numHours = int(totTimeSec) // 3600
    numMins = int(totTimeSec % 3600) // 60
    numSecs = (totTimeSec % 3600) % 60
    if (continuation == 1):
        f.write("ADDED running time: "+str(numDays)+ " days "+str(numHours)+ " hours "+str(numMins)+" minutes "+str(numSecs)+" seconds.\n")
    else:
        f.write("Running time: "+str(numDays)+ " days "+str(numHours)+ " hours "+str(numMins)+" minutes "+str(numSecs)+" seconds.\n")

    # close the file
    f.close()



# Call main()
main()
