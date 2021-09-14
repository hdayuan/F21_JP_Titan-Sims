import rebound
import numpy as np
import matplotlib.pyplot as plt

G = 6.67e-11
M_SUN = 1.989e30

sim = rebound.Simulation()
sim.units = ('AU', 'days', 'Msun')

# add Saturn: 0.0002857 solar masses
mSaturn = 0.0002857
sim.add(m=mSaturn)

# add Titan: 6.76319759 * 10^-8 solar masses, a AUs semi-major axis, and velocity
# based on a
# Note: Assuming perfectly circular orbit to begin with, so a is constant
a = 0.008167696467
v = np.sqrt(mSaturn * G * M_SUN / (a * 1.496e11))
sim.add(m=0.0000000676319759, x=a, vy=v)

sim.move_to_com()
fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True, periastron=True)

sim.integrator = "whfast"
sim.dt = 1e-3
sim.integrate(6.28318530717959, exact_finish_time=0)   # 6.28318530717959 is 2*pi

particles = sim.particles
torb = 2.*np.pi
Noutputs = 100
times = np.linspace(torb, 2.*torb, Noutputs)
x = np.zeros(Noutputs)
y = np.zeros(Noutputs)
for i,time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
plt.plot(x, y)





"""# We can add Jupiter and four of its moons by name, since REBOUND is linked to the HORIZONS database.
labels = ["Jupiter", "Io", "Europa","Ganymede","Callisto"]
sim.add(labels)

os = sim.calculate_orbits()
print("n_i (in rad/days) = %6.3f, %6.3f, %6.3f" % (os[0].n,os[1].n,os[2].n))
print("P_i (in days)     = %6.3f, %6.3f, %6.3f" % (os[0].P,os[1].P,os[2].P))

sim.move_to_com()
fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=True, periastron=True)

sim.integrator = "whfast"
sim.dt = 0.05 * os[0].P  # 5% of Io's period
Nout = 100000            # number of points to display
tmax = 800*365.25         # let the simulation run for 80 years
Nmoons = 4"""
