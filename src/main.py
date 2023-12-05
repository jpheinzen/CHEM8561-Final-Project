import numpy as np
# import matplotlib.pyplot as plt
from functions import *
import time
import sys

file = open("/project/heinz194/private/classes/CHEM_8561/Project/U.txt", "w")


# 3d plot, how to use VScode after plotting


# rng = np.random.default_rng()
bg = np.random.MT19937(1)
rng = np.random.Generator(bg)


# Initializing points in box

nPts = 108
# box = 2e-8
# dp = 5e-9
# V = box**3
# dv = V/2**3
# print(V,dv,(V+dv)**(1/3),box)
T = 320
P = 1e5
P = 7538580
P = 1.005e+8
k =  1.380649e-23        # J/K.
beta = 1/(k*T)
nAccepts = 0
ntot = 0
nVAccepts = 0
nVtot = 0

V = nPts*k*T/P
box = V**(1/3)
dp = box
dv = V/2**3

N = 51
debug = False
# debug = True

# stuff to write to file
Usave = np.zeros(N*nPts)
iUsave = 0
file.write('%i\t%i\n' % (N,nPts))


tic = time.perf_counter()


# Generate configuration
pts = generateConfig(nPts, box, rng)

# plotParticles(pts,box,box)


# Calculate energy Uold
(LJ12old,LJ6old) = getPotentials(pts, nPts, box)
Uold = np.sum(LJ12old+LJ6old)
# Uold = getPotentialOld(pts, nPts, box)

# print((Uold-UU)/UU)
# sys.exit()

for n in range(N):
    print()
    print('----------------------------- iter',n,'-----------------------------')

    for i in range(nPts):
        # print(i,'iter')
        if debug:
            # FIGURE OUT WHY THIS DOENS'T WORK???
            # (LJ12old2,LJ6old2) = getPotentials(pts, nPts, box)
            # Uold2 = np.sum(LJ12old2+LJ6old2)
            Uold2 = getPotentialOld(pts, nPts, box)
            if abs((Uold2-Uold)/Uold) > 1e-3:
                print('Uold2',Uold,'Uold',Uold2,'diff',Uold-Uold2)
                print('pt1',pts[1,1])
                sys.exit()

        # displace particle coords

        # oldPt = pts[:,i:i+1]
        shift = dp*(rng.random([3,1])-0.5)
        pts[:,i:i+1] += shift
        pts[:,i:i+1] %= box
        # print('shift',shift)
        # print('oldPt',oldPt)

        # Calculate Un
        (LJ12newsmall,LJ6newsmall) = getPotential(pts, i, box)

        I = getInds(i,nPts)
        dU = np.sum(LJ12newsmall+LJ6newsmall) - np.sum(LJ12old[I]+LJ6old[I])
        Unew = getPotentialOld(pts, nPts, box)
        # Unew = np.sum(LJ12new+LJ6new)

        # print('old',Uold,'\tnew',Unew,'\tdiff',Unew-Uold)
        print("old: %#.5g \t new: %#.5g \t diff:% #.5g \t box: %#.5g" % (Uold, Unew, Unew-Uold,box))
        # print("old: %#.5g \t new: %#.5g \t diff:% #.5g \t box: %#.5g" % (Uold, Uold+dU, dU,box))

        print(Unew-Uold,dU,Unew-Uold-dU)

        # Accept/reject move
        if (Unew-Uold) < 0 or acceptMove(Unew-Uold, beta, rng):
            nAccepts += 1
            LJ12old[I] = LJ12newsmall
            LJ6old[I] = LJ6newsmall
            Uold = Unew
            # print('passed')
        else:
            pts[:,i:i+1] -= shift
            # print('pt',pts[:,i:i+1],'vs old', oldPt)
            # print('failed')

        ntot += 1
        Usave[iUsave] = Uold
        file.write("%i\t%g\n" % (iUsave, Usave[iUsave]))
        iUsave +=1
     
    # change volume
    Vnew = V + dv*(rng.random()-0.5)
    boxNew = Vnew**(1/3)
    coordFact = boxNew/box
    pts *= coordFact
    (LJ12new,LJ6new) = getPotentials(pts, nPts, boxNew)
    Unew = np.sum(LJ12new+LJ6new)

    dH = Unew-Uold + P*dv - k*T*nPts*np.log(Vnew/V) # type: ignore
    print('\n\n\n\n')
    # print('dV',Vnew-V,'boxdiff',boxNew-box,)
    print("boxNew: %#.5g \t box: %#.5g \t dH: % #.5g" % (boxNew, box, dH))

    # Accept/Reject Volume Change
    if (dH) < 0 or acceptMove(dH, beta, rng):
        nVAccepts += 1
        V = Vnew
        Uold = Unew
        box = boxNew
        LJ12old = LJ12new
        LJ6old = LJ6new
        print('passed')
    else:
        pts /= coordFact
        # pass
        # pts[:,i:i+1] -= shift
        print('failed')
    nVtot += 1

    if n % 10  == 0:
        percPass = nAccepts/ntot
        print(percPass*100, '% passed','\td', dp)
        if percPass > 0.5:
            dp *= 1.05
        elif percPass < 0.5:
            dp *= 0.95
        ntot = 0
        nAccepts = 0


        # plotParticles(pts,box,box)
        

        # plotParticles(ptsNew,boxNew,boxNew)
        if n % 10 == 0:     
          
            # Try to keep the percentage of Volume moves that pass to 50%
            percVPass = nVAccepts/nVtot
            print(percVPass*100, '% passed','\tdv', dv)
            if percVPass > 0.5:
                dv *= 1.05
            elif percVPass < 0.5:
                dv *= 0.95
            
            nVAccepts = 0
            nVtot = 0

            # time.sleep(5)

            
        print('\n\n\n\n')
       

    # if n % N/1000 == 0:
    #     plotParticles(pts)


toc = time.perf_counter()
totalTime = toc-tic
print(f"Elapsed time: {totalTime:0.4f} seconds")

# plt.plot(np.arange(N*nPts),Usave,'o')

# plotParticles(pts,box,box)

Uold2 = getPotentialOld(pts, nPts, box)
print(abs((Uold2-Uold)/Uold),Uold2)

file.close()

# potentials = getPotentials((np) sqDistance, bool isCO2)
# accept = acceptMove(double dPotential)


# old: -8.2495e-19         new: -4.5717e-19        diff: 3.6777e-19        box: 2.1008e-09

# boxNew: 2.0834e-09       box: 2.1008e-09         dH:  1.3635e-19
# failed
# 0.6481481481481481 % passed     d 1.3005139389942589e-09
# 30.0 % passed   dv 5.90509347176418e-28

# Elapsed time: 15.4381 seconds


# 1. Generate configuration
# 2. Calculate energy Uold
# 3. displacement of particle coords
# 4. Calculate Un
# 5. accept/reject move
# 	if accepted, keep
# 	if rejected, undo move
# 6. every X iterations, change volume

# repeat at 3
