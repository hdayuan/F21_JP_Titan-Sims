import sys
import numpy as np


timescale = 810766627.6731131 # in years
iaTitanRS = float(sys.argv[1])
target_runtime = float(sys.argv[2])
h_or_d = sys.argv[3]

if h_or_d == 'd':
    # target_runtime in days, convert to hours
    target_runtime = target_runtime*24.

# 87 hours 6 minutes 29.20167827606201 seconds runtime for 192696243.6859105 year sim
# number of years sim can run in 1 hour (years/hour)
runT_TO_simT = 192696243.6859105/(87.+(6./60.)+(29.20167827606201/3600.))

faTitanRS = iaTitanRS*np.exp(target_runtime*runT_TO_simT/timescale)

print(faTitanRS)
