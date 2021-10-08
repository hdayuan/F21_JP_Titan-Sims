import rebound
import numpy as np
import matplotlib.pyplot as plt


sim = rebound.Simulation()
sim.units = ('AU', 'days', 'Msun')

rebound.data.add_solar_system(sim)

"""sim.add(m=1, x=0, y=0, vx=0, vy=0) # add sun
sim.add(m=1e-3, a=2., e=0.1) # add planet"""

sim.integrate(10)
# move to helio centric frame of reference
sim.move_to_hel()

os = sim.calculate_orbits()
print("n_i (in rad/days) = %6.3f, %6.3f, %6.3f" % (os[0].n,os[1].n,os[2].n))
print("P_i (in days)     = %6.3f, %6.3f, %6.3f" % (os[0].P,os[1].P,os[2].P))

