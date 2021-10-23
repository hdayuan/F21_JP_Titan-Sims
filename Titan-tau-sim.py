"""Description:
Takes 4 command line arguments: number of samples, initial and 
final semi-major axes of titan in Saturn radii, and factor to divide
timescale by
Prints these four parameters in the order above,
then calls main() (see below for description of main()). Finally, prints out
total time for the program to run."""

import sys
import rebound
import reboundx
import numpy as np
import time


"""4 parameters: numSamples is number of samples over integration, 
ia_titan and fa_titan are initial and final semi-major axes of Titan 
in Saturn radii respectively, file is the file to which output
is written

Integrates the system from a = ia_titan to a = fa_titan and prints 
numSamples data points, each on one line in the following format: 

'semi-major axis in Saturn radii'[tab]'eccentricity'[tab]'longitude of pericenter'
[tab]'mean anomaly of Sun's "orbit" around Saturn'[tab]'current time in simulation'
"""
def main(numSamples, ia_titanRS, fa_titanRS, file, factor):

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

    ia_titanAU = ia_titanRS * rSat  # starting semi-major axis of Titan

    # Calculate the expected semi-major axis of Titan at which
    # evection resonance should occur given i = e = 0
    # in units of saturn radii
    exp_aResRS = np.power((9./4.)*j2Sat*j2Sat*mSat*(aSat/rSat)**3, (1./7.))

    # calculate period of Titan in years at evection resonance
    # distance ** Check accuracy **
    mm = np.sqrt(G*mSat*mSun / ((exp_aResRS*rSat*AU_TO_M)**3))  # in rad / sec
    tauTitan = 2*np.pi*(1/mm)/(3600*24*365.25) # in years

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'Msun')
    sim.integrator = "whfast"
    sim.dt = (1./20.) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat, hash = "Saturn")

    # add Titan
    sim.add(m=mTitan, a=ia_titanAU, e=eTitan, hash = "Titan")

    # add sun (with semi-major axis and eccentricity of Saturn)
    sim.add(m=1., a=aSat, e=eSat, hash = "Sun")

    # Initiate reboundx
    rebx = reboundx.Extras(sim)

    # Calculate timescale for exponential migration of Titan's semi-major axis
    ageSat = 4.503e9 # age of saturn in yrs
    timescale = (1./factor) * 3. * ((exp_aResRS*rSat/curr_aTitan)**3) * ageSat

    # calculate total time to integrate based on ia_titan, fa_titan, tau, and
    # given exponential migration of form a = a_0 * e^(t/tau)
    totSimTime = timescale * np.log(fa_titanRS / ia_titanRS) 

    # Add migration force for Titan's outward migration (a = a0e^(t/tau)
    mof = rebx.load_force("modify_orbits_forces")
    rebx.add_force(mof)
    sim.particles["Titan"].params["tau_a"] = timescale

    # add Saturn's J2
    gh = rebx.load_force("gravitational_harmonics")
    rebx.add_force(gh)
    sim.particles["Saturn"].params["J2"] = j2Sat
    sim.particles["Saturn"].params["R_eq"] = rSat

    plotDT = totSimTime/numSamples

    # Integrate
    for i in range(numSamples):
        sim.integrate(i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(sim.particles["Titan"].a/rSat)+"\t"+str(sim.particles["Titan"].e)+"\t")
        file.write(str(sim.particles["Titan"].pomega)+"\t"+str(sim.particles["Sun"].l))
        file.write("\t"+str(sim.t)+"\n")

    # sim.cite()



# start timer
start_time = time.time()

# Parse command-line arguments
numSamples = int(sys.argv[1])
ia_titan = float(sys.argv[2])
fa_titan = float(sys.argv[3])
div_factor = float(sys.argv[4])

# open file
f = open("v2-div"+str(div_factor)+"-output-"+str(numSamples)+"s-"+str(ia_titan)+"to"+str(fa_titan)+"rs.txt", "a")

# Write parameters of simulation
f.write(str(numSamples)+"\n")
f.write(str(ia_titan)+"\n")
f.write(str(fa_titan)+"\n")
f.write(str(div_factor)+"\n")

# call main
main(numSamples, ia_titan, fa_titan, f, div_factor)

# Write running time
totTimeSec = time.time() - start_time
numHours = int(totTimeSec) // 3600
numMins = int(totTimeSec % 3600) // 60
numSecs = (totTimeSec % 3600) % 60
f.write("Running time: "+str(numHours)+ " hours "+str(numMins)+" minutes "+str(numSecs)+" seconds.\n")

# close the file
f.close()

