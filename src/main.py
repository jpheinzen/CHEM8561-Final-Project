import numpy as np
import matplotlib.pyplot as plt
from functions import *

# 3d plot, how to use VScode after plotting


# rng = np.random.default_rng()
bg = np.random.MT19937(1)
rng = np.random.Generator(bg)


# Initializing points in box

nPts = 108
box = 2e-8
# box = 5e-7
dp = box/5
# dp = 5e-9
V = box**3
dv = V/2**3
# print(V,dv,(V+dv)**(1/3),box)
T = 500
P = 1e5
k =  1.380649e-23        # J/K.
beta = 1/(k*T)
nAccepts = 0
ntot = 0
nVAccepts = 0
nVtot = 0
N = 1000

V2 = nPts*k*T/P
print(V,V2,V2**(1/3),box)


# Generate configuration
pts = generateConfig(nPts, box, rng)

plotParticles(pts,box,box)


# Calculate energy Uold
Uold = getPotential(pts, nPts, box)
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
        Unew = getPotential(pts, nPts, box)
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

    if n % 100 == 0:
        # plotParticles(pts,box,box)
        # change volume
        ptsNew = pts
        Vnew = V + dv*(rng.random()-0.5)
        boxNew = Vnew**(1/3)
        coordFact = boxNew/box
        ptsNew *= coordFact
        Unew = getPotential(ptsNew, nPts, boxNew)

        dH = Unew-Uold + P*dv - k*T*nPts*np.log(Vnew/V)
        print('\n\n\n\n')
        # print('dV',Vnew-V,'boxdiff',boxNew-box,)
        print("boxNew: %#.5g \t box: %#.5g \t dH: % #.5g" % (boxNew, box, dH))

        # plotParticles(ptsNew,boxNew,boxNew)

        # Accept/Reject Volume Change
        if (dH) < 0 or acceptMove(dH, beta, rng):
            nVAccepts += 1
            pts = ptsNew
            V = Vnew
            Uold = Unew
            box = boxNew
            print('passed')
        else:
            pass
            # pts[:,i:i+1] -= shift
            print('failed')
        nVtot += 1

        if n % 100 == 0:
            # Try to keep the percentage of Volume moves that pass to 50%
            percVPass = nVAccepts/nVtot
            print(percVPass*100, '% passed','\td', dv)
            if percVPass > 0.5:
                dv *= 1.05
            elif percVPass < 0.5:
                dv *= 0.95
        print('\n\n\n\n')

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
