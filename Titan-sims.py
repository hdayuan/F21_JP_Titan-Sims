"""Description:
Takes 4 command line arguments: number of samples, total simulation time
in years, initial semi-major axis of titan in Saturn radii, and timescale
for Titan's migration. Prints these four parameters in the order above,
then calls main() (see below for description of main()). Finally, prints out
total time for the program to run."""

import sys
sys.path.append('/u/hdyuan/bin')

""" TO-DO:
- Plot expected precession period
"""

import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt
import time


"""4 parameters: numSamples is number of samples over integration, 
totalSimTime is total time of simulation in years, ia_titan is initial
semi-major axis of Titan in Saturn radii, timescale is the timescale (tau)
for the migration of Titan's semi-major axis

Integrates the system for totalSimTime years and prints numSamples
data points, each in the following order: 
    semi-major axis in Saturn radii
    eccentricity
    current time

If timescale is 0, calculates correct timescale
"""
def main(numSamples, totalSimTime, ia_titan, timescale):

    # Constants
    G = 6.67e-11 # G in SI units ****
    mSun = 1.9891e30 # mass of sun in kg ****
    AU_TO_M = 1.496e+11 # meters in one AU

    # Check accuracy of all the following constants
    mSat = 0.0002857 # mass of saturn in solar masses
    aSat = 9.5549 # semi-major axis of Saturn's orbit in AU
    eSat = 0.0565 # eccentricity of Saturn
    rSat = 0.00038926024 # radius of Saturn in AU
    j2Sat = 16298e-6 # J2 of Saturn (Murray and Dermott p 531) 
    curr_aTitan = 0.008167696467 # modern-day semi-major axis of titan in AU
    mTitan = 0.0000000676319759 # mass of titan in solar masses
    eTitan = 0.001 # eccentricity of Titan's orbit ***** Modify as needed *****

    ini_aTitan = ia_titan * rSat  # starting semi-major axis of Titan

    # Calculate the expected semi-major axis of Titan at which
    # evection resonance should occur given i = e = 0
    # in units of saturn radii
    exp_aRes = np.power((9/4)*j2Sat*j2Sat*mSat*(aSat/rSat)**3, (1/7))

    # calculate period of Titan in years at evection resonance
    # distance ** Check accuracy **
    mm = np.sqrt(G*mSat*mSun / ((exp_aRes*rSat*AU_TO_M)**3))  # in rad / sec
    tauTitan = 2*np.pi*(1/mm)/(3600*24*365.25) # in years

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'Msun')
    sim.integrator = "whfast"
    sim.dt = (1/20) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat)

    # add Titan
    sim.add(m=mTitan, a=ini_aTitan, e=eTitan)

    # add sun (with semi-major axis and eccentricity of Saturn)
    sim.add(m=1, a=aSat, e=eSat)

    # Initiate reboundx
    rebx = reboundx.Extras(sim)

    # Calculate timescale for exponential migration of Titan's semi-major axis
    if (timescale == 0):
        ageSat = 4.503e9 # age of saturn in yrs
        tau = 3 * ((exp_aRes*rSat/curr_aTitan)**3) * ageSat
    else:
        tau = timescale

    # Add migration force for Titan's outward migration (a = a0e^(t/tau)
    mof = rebx.load_force("modify_orbits_forces")
    rebx.add_force(mof)
    sim.particles[1].params["tau_a"] = tau

    # add Saturn's J2
    gh = rebx.load_force("gravitational_harmonics")
    rebx.add_force(gh)
    sim.particles[0].params["J2"] = j2Sat
    sim.particles[0].params["R_eq"] = rSat

    plotDT = totalSimTime/numSamples

    # Integrate
    for i in range(0,numSamples):
        sim.integrate(i * plotDT)
        print(sim.particles[1].a / rSat)
        print(sim.particles[1].e)
        # print(sim.particles[1].pomega)
        print(sim.t)

    # sim.move_to_hel()




start_time = time.time()

numSamples = int(sys.argv[1])
totalSimTime = int(sys.argv[2])
ia_titan = float(sys.argv[3])
timescale = int(sys.argv[4])

# Print parameters of simulation
print(numSamples)
print(totalSimTime)
print(ia_titan)
print(timescale)

main(numSamples, totalSimTime, ia_titan, timescale)
print(str(time.time() - start_time) + " seconds")
