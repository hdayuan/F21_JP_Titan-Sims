"""Description:

Questions: mass (M) and mean motion (n) in T to Q conversion --> check plot timescale too

To run:

python3 Titan-sims.py [number of samples] [initial a in Saturn radii]
[target final a in Saturn radii]

Takes 3 command line arguments: number of samples, initial and
final semi-major axes of titan in Saturn radii
Prints these four parameters in the order above,
then calls main() (see below for description of main()). Finally, prints out
total time for the program to run."""

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

"""Calculates tidal time-lag T (in years) for a body with mean motion n (radians/sec), 
tidal Q factor Q, radius r (AU), and mass m (solar masses)"""
def time_lag(n, Q, r, m):
    return 2*n*Q*(r*AU_TO_M)**3 / (G*m*M_SUN*YR_TO_SEC)

"""Calculates the mean motion (in radians/sec) of a body orbiting around a central body
of mass mCentral (solar masses) at semi-major axis a (AU) """
def mean_motion(mCentral, a):
    return np.sqrt(G*mCentral*M_SUN / ((a*AU_TO_M)**3))

"""4 parameters: numSamples is number of samples over integration,
ia_titan and fa_titan are initial and final semi-major axes of Titan
in Saturn radii respectively, file is the file to which output
is written

Integrates the system from a = ia_titan to a = fa_titan and prints
numSamples data points, each on one line in the following format:

'semi-major axis in Saturn radii'[tab]'eccentricity'[tab]'longitude of pericenter'
[tab]'mean anomaly of Sun's "orbit" around Saturn'[tab]'current time in simulation'
"""
def main(numSamples, ia_titanRS, fa_titanRS, file):

    # Check accuracy of all the following constants
    mSat = 0.0002857 # mass of saturn in solar masses
    aSat = 9.5549 # semi-major axis of Saturn's orbit in AU
    eSat = 0.0565 # eccentricity of Saturn
    rSat = 0.00038926024 # radius of Saturn in AU
    oSat = 26.7 * np.pi / 180. # obliquity of Saturn in radians
    j2Sat = 16298e-6 # J2 of Saturn (Murray and Dermott p 531)
    # calculated later # QSat = 5000. # Tidal Q factor of Saturn (Lainey et al.)
    k2Sat = 0.39 # Love number of Saturn *** CHECK THIS ***
    curr_aTitan = 0.008167696467 # modern-day semi-major axis of titan in AU
    mTitan = 0.0000000676319759 # mass of titan in solar masses
    eTitan = 0.001 # initial eccentricity of Titan's orbit ***** Modify as needed *****
    rTitan = 0.04421567543 * rSat # radius of Titan in AU
    k2Titan = 0.15 # Love number of Titan *** CHECK THIS ***
    QTitan = 100. # Tidal Q factor of Titan

    ia_titanAU = ia_titanRS * rSat  # starting semi-major axis of Titan

    # Calculate the expected semi-major axis of Titan at which
    # evection resonance should occur given i = e = 0
    # in units of saturn radii
    exp_aResRS = np.power((9./4.)*j2Sat*j2Sat*mSat*(aSat/rSat)**3, (1./7.))

    # calculate period of Titan in years at evection resonance
    # distance ** Check accuracy **
    nTitan = mean_motion(mSat, exp_aResRS*rSat) # in rad / sec
    tauTitan = 2*np.pi*(1/nTitan)/YR_TO_SEC # in years

    # Initialize rebound simulation
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'MSun')
    sim.integrator = "whfast"
    sim.dt = (1./20.) * tauTitan # time step = 1/20 * shortest orbital period

    # add Saturn
    sim.add(m=mSat, hash = "Saturn")

    # add Titan
    sim.add(m=mTitan, a=ia_titanAU, e=eTitan, hash = "Titan")

    # add Sun (with semi-major axis and eccentricity of Saturn, inclination = obliquity of Saturn)
    sim.add(m=1., a=aSat, e=eSat, hash = "Sun")

    # Calculate timescale for exponential migration of Titan's semi-major axis
    ageSat = 4.503e9 # age of saturn in yrs
    timescale = 3. * ((exp_aResRS*rSat/curr_aTitan)**3) * ageSat

    # calculate total time to integrate based on ia_titan, fa_titan, tau, and
    # given exponential migration of form a = a_0 * e^(t/tau)
    totSimTime = timescale * np.log(fa_titanRS / ia_titanRS)

    # Initiate reboundx
    rebx = reboundx.Extras(sim)
    ps = sim.particles

    # Add migration force for Titan's outward migration (a = a0e^(t/tau)
    # mof = rebx.load_force("modify_orbits_forces")
    # rebx.add_force(mof)
    #ps["Titan"].params["tau_a"] = timescale

    # add Saturn's J2
    gh = rebx.load_force("gravitational_harmonics")
    rebx.add_force(gh)
    ps["Saturn"].params["J2"] = j2Sat
    ps["Saturn"].params["R_eq"] = rSat

    # add tidal forces
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    # Tidal forces of Titan
    titanT = time_lag(nTitan, QTitan, rTitan, mTitan) # in years
    ps["Titan"].r = rTitan # AU
    ps["Titan"].params["tctl_k2"] = k2Titan
    ps["Titan"].params["tctl_tau"] = titanT
    ps["Titan"].params["Omega"] = nTitan*YR_TO_SEC # in radians per year

    # Tidal forces of Saturn
    nSat = mean_motion(1., aSat)  # in rad / sec
    QSat = 3*k2Sat*(mTitan/mSat)*nTitan*timescale*(1.0/exp_aResRS)**5. # Q for Saturn 

    satT = time_lag(nSat, QSat, rSat, mSat) # in years
    ps["Saturn"].r = rSat # AU
    ps["Saturn"].params["tctl_k2"] = k2Sat
    ps["Saturn"].params["tctl_tau"] = satT
    omega = 2000000*np.pi*YR_TO_SEC/(10.656*3600)
    ps["Saturn"].params["Omega"] = omega # in rad per year

    plotDT = totSimTime/numSamples

    # Integrate
    for i in range(numSamples):
        sim.integrate(i * plotDT)
        # move to Saturn's frame of reference
        sim.move_to_hel()
        # write output
        file.write(str(ps["Titan"].a/rSat)+"\t"+str(ps["Titan"].e)+"\t")
        file.write(str(ps["Titan"].pomega)+"\t"+str(ps["Sun"].l)+"\t")
        file.write(str(ps["Titan"].inc)+"\t")
        file.write(str(sim.t)+"\n")

    # sim.cite()
    return totSimTime

# start timer
start_time = time.time()

# Parse command-line arguments
numSamples = int(sys.argv[1])
ia_titan = float(sys.argv[2])
fa_titan = float(sys.argv[3])

# open file
f = open("v4out-tides-"+str(numSamples)+"s-"+str(ia_titan)+"to"+str(fa_titan)+"rs.txt", "a")

# Write parameters of simulation
f.write(str(numSamples)+"\n")
f.write(str(ia_titan)+"\n")
f.write(str(fa_titan)+"\n")

# call main
totSimTime = main(numSamples, ia_titan, fa_titan, f)

# Write total sim time
f.write("Total Integration Time: "+str(totSimTime)+" years.\n")

# Write running time
totTimeSec = time.time() - start_time
numHours = int(totTimeSec) // 3600
numMins = int(totTimeSec % 3600) // 60
numSecs = (totTimeSec % 3600) % 60
f.write("Running time: "+str(numHours)+ " hours "+str(numMins)+" minutes "+str(numSecs)+" seconds.\n")

# close the file
f.close()
