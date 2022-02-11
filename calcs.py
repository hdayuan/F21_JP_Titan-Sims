
AU_TO_M = 1.496e+11 # meters in one AU
rSat = 58232503./AU_TO_M # mean physical radius
aIap = 3561000000.
aTitan = 0.008167696467 # modern-day semi-major axis
iaTitan = 8.2 * rSat  # starting semi-major axis of Titan
t_start = (iaTitan/aTitan)**3
iaIap = aIap * (t_start)**(1./3.) # initial semi-major axis of Iapetus
print(iaIap)

