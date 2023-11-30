import numpy as np
import matplotlib.pyplot as plt
from functions import *

# 3d plot, how to use VScode after plotting


# rng = np.random.default_rng()
bg = np.random.MT19937(1)
rng = np.random.Generator(bg)


# Initializing points in box

nPts = 108
box = 5e-8
# box = 5e-7
dp = box/5
# dp = 5e-9
V = box**3
dv = V/5**3
hbox = box/2
T = 500
P = 1e5
k =  1.380649e-23        # J/K.
beta = 1/(k*T)
nAccepts = 0
ntot = 0
N = 10

# Generate configuration
pts = generateConfig(nPts, box, rng)

plotParticles(pts,box,box)


# Calculate energy Uold
Uold = getPotential(pts, nPts, box, hbox)
print()

for n in range(N):
    print()
    print('----------------------------- iter',n,'-----------------------------')

    for i in range(nPts):
        # print(i,'iter')

        # displace particle coords

        # oldPt = pts[:,i:i+1]
        shift = dp*(rng.random([3,1])-0.5)
        pts[:,i:i+1] += shift
        pts[:,i:i+1] %= box
        # print('shift',shift)
        # print('oldPt',oldPt)

        # Calculate Un
        Unew = getPotential(pts, nPts, box, hbox)
        # print('old',Uold,'\tnew',Unew,'\tdiff',Unew-Uold)
        print("old: %#.5g \t new: %#.5g \t diff: % #.5g" % (Uold, Unew, Unew-Uold))

        # Accept/reject move
        if (Unew-Uold) < 0 or acceptMove(Unew-Uold, beta, rng):
            nAccepts += 1
            Uold = Unew
            # print('passed')
        else:
            pts[:,i:i+1] -= shift
            # print('pt',pts[:,i:i+1],'vs old', oldPt)
            # print('failed')

        ntot += 1

    # Try to keep the percentage of moves that pass to 50%
    percPass = nAccepts/ntot
    print(percPass*100, '% passed','\td', dp)
    if percPass > 0.5:
        dp *= 1.05
    elif percPass < 0.5:
        dp *= 0.95

    # if n % N/1000 == 0:
    #     plotParticles(pts)




plotParticles(pts,box,box)




# potentials = getPotentials((np) sqDistance, bool isCO2)
# accept = acceptMove(double dPotential)



# 1. Generate configuration
# 2. Calculate energy Uold
# 3. displacement of particle coords
# 4. Calculate Un
# 5. accept/reject move
# 	if accepted, keep
# 	if rejected, undo move
# 6. every X iterations, change volume

# repeat at 3
