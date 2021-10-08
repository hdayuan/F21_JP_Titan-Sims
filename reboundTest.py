import rebound
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.units = ('AU', 'days', 'Msun')

titanT = 16 # period of titan in days

# add Saturn: 0.0002857 solar masses
mSaturn = 0.0002857
sim.add(m=mSaturn)

# add Titan: 6.76319759 * 10^-8 solar masses, a AUs semi-major axis, and velocity
# based on a
# Note: Assuming approximately circular orbit (e = 0.001) to begin with, so a is constant
# titan = rebound.Particle(G=sim.G, primary=saturn, m=0.0000000676319759, a=0.008167696467, e=0.1)
aTitan = 0.008167696467
mTitan = 0.0000000676319759
eTitan = 0.1
sim.add(m=mTitan, a=aTitan, e=eTitan)

# add sun
aSat = 9.5549
eSat = 0.0565
sim.add(m=1, a=aSat, e=eSat)

sim.integrator = "whfast"
sim.dt = (1/20) * titanT # time step = 1/20 * shortest orbital period
sim.integrate(20 * 365)   # simulate for 80 years

sim.move_to_hel()

fig = rebound.OrbitPlot(sim, unitlabel="[AU]", color=False, periastron=True)

ax = plt.subplot(111)
ax.set_xlim([-10,10])
ax.set_ylim([-10,10])

plt.show()

"""particles = sim.particles
torb = 2.*np.pi
Noutputs = 10
times = np.linspace(torb, 2.*torb, Noutputs)
x = np.zeros(Noutputs)
y = np.zeros(Noutputs)
for i,time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y"""
