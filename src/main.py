import numpy as np
# import matplotlib.pyplot as plt
from functions import *
import time
import sys

# rng = np.random.default_rng()
bg = np.random.MT19937(1)
rng = np.random.Generator(bg)

atm2Pa = 101325

# Initializing points in box
nPts = 108
# nPts = 5
T = 320
P = 74.4*atm2Pa
P = 992*atm2Pa
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
N = 1_000_000+5001
N = 5001
debug = False
# debug = True
writeUTF = False
writeUTF = True
writeptsTF = True
writeEvery = 100
writeToScreen = False
writeToScreen = True

if writeUTF: # energy to write to file
    file1 = open("U.txt", "w")
    # Usave = np.zeros(N*nPts)
    # iUsave = 0
    iTot = 0
    file1.write('%i\t%i\t%i\n' % (N,nPts,writeEvery))
if writeptsTF:  # points to write to file
    file2 = open("pts.txt", "w")
    file2.write('%i\t%i\t%i\n' % (N,nPts,writeEvery))



tic = time.perf_counter()

# Generate configuration
pts = generateConfig(nPts, box, rng)

# plotParticles(pts,box,box)


# Calculate energy Uold
(LJ12old,LJ6old) = getPotentials(pts, nPts, box)
Uold = addLJArr(LJ12old,LJ6old)
# UU = getPotentialOld(pts, nPts, box)

# print((Uold-UU)/UU,Uold,UU)
# sys.exit()

for n in range(N):
    if writeToScreen or n % int(min(N/10,1000)) == 0:
        print()
        print('----------------------------- iter',n,'-----------------------------')

    for i in range(nPts):
        # print(i,'iter')
        if debug:
            check(Uold,getPotentialOld(pts, nPts, box),'1',tol = 5e-12)

        # displace particle coords

        # oldPt = pts[:,i:i+1]
        shift = dp*(rng.random([3,1])-0.5)
        # oldPt = np.copy(pts[:,i:i+1])
        pts[:,i:i+1] += shift
        pts[:,i:i+1] %= box
        # print('shift',shift)
        # print('oldPt',oldPt)

        # Calculate Un
        (LJ12newsmall,LJ6newsmall) = getPotential(pts, i, box)

        I = getInds(i,nPts)
        dU = addLJArr(LJ12newsmall,LJ6newsmall) - addLJArr(LJ12old[I],LJ6old[I])
        
        if debug:
            check(Uold+dU,getPotentialOld(pts, nPts, box),i,exit=False,tol=1e-12)

        # print("old: %#.5g \t new: %#.5g \t diff:% #.5g \t box: %#.5g" % (Uold, Unew, Unew-Uold,box))
        if writeToScreen:
            print("old: %#.5g \t new: %#.5g \t diff:% #.5g \t box: %#.5g" % (Uold, Uold+dU, dU,box))

        # Accept/reject move
        if (dU) < 0 or acceptMove(dU, beta, rng):
            nAccepts += 1
            
            Uold += dU
            # Uold = Unew

            LJ12old[I] = LJ12newsmall
            LJ6old[I] = LJ6newsmall
            # print('passed')
        else:
            pts[:,i:i+1] -= shift
            pts[:,i:i+1] %= box
            # print('failed')

        ntot += 1

        if writeUTF and ((iTot % (nPts*writeEvery) == 0) or (iTot < writeEvery*nPts)):
        # if writeUTF and ((iTot % (nPts*writeEvery) == 0)):
            # Usave[iUsave] = Uold
            file1.write("%i\t%g\n" % (iTot, Uold))
            # iUsave +=1
        iTot += 1


    if writeptsTF and (n % writeEvery == 0):
        file2.write('%i\t%g\t%g\t%g\n' % (n,V,dp,dv))
        writePts(file2,pts,nPts)

    # change volume
    Vnew = V + dv*(rng.random()-0.5)
    boxNew = Vnew**(1/3)
    volFact = V/Vnew
    LJ12new = LJ12old*volFact**4
    LJ6new = LJ6old*volFact**2
    coordFact = boxNew/box
    # pts *= coordFact
    # (b,c) = getPotentials(pts, nPts, boxNew)
    # Unew3 = 4*eps*np.sum(b-c)
    # check(LJ12new,b,'12')
    # check(LJ6new,c,'6')

    Unew = addLJArr(LJ12new,LJ6new)

    if debug:
        check(Unew,getPotentialOld(pts*coordFact, nPts, boxNew),'3',tol=5e-15)
    

    dH = Unew-Uold + P*dv - k*T*nPts*np.log(Vnew/V) # type: ignore
    if writeToScreen:
        print('\n\n\n\n')
        print("boxNew: %#.5g \t box: %#.5g \t dH: % #.5g" % (boxNew, box, dH))

    # Accept/Reject Volume Change
    if (dH) < 0 or acceptMove(dH, beta, rng):
        nVAccepts += 1
        pts *= coordFact
        V = Vnew
        Uold = Unew
        box = boxNew
        LJ12old = LJ12new
        LJ6old = LJ6new
        
        if writeToScreen:
            print('passed')
    else:
        # pts /= coordFact
        if writeToScreen:
            print('failed')
    nVtot += 1

    if n % 10  == 0:
        # Try to keep the percentage of particle moves that pass to 50%
        percPass = nAccepts/ntot
        if percPass > 0.5:
            dp = min(box,dp*1.05)
        elif percPass < 0.5:
            dp *= 0.95
        ntot = 0
        nAccepts = 0


        # Try to keep the percentage of Volume moves that pass to 50%
        percVPass = nVAccepts/nVtot
        if percVPass > 0.5:
            dv *= 1.05
        elif percVPass < 0.5:
            dv *= 0.95
        
        nVAccepts = 0
        nVtot = 0


        if writeToScreen:
            print(percPass*100, '% passed','\td', dp)
            print(percVPass*100, '% passed','\tdv', dv)
            print('\n\n\n\n')
       

    # if n % N/1000 == 0:
    #     plotParticles(pts)

# plt.plot(np.arange(N*nPts),Usave,'o')

# plotParticles(pts,box,box)
toc = time.perf_counter()
totalTime = toc-tic
print(f"Elapsed time: {totalTime:0.4f} seconds")
totalTime /= 3600
print(f"Elapsed time: {totalTime:0.4f} hours")

print(percPass*100, '% passed','\td', dp)
print(percVPass*100, '% passed','\tdv', dv)

Uold2 = getPotentialOld(pts, nPts, box)
print('Uold vs Uold2 relative:',abs((Uold2-Uold)/Uold),'\tUold2',Uold2)
print('done')

if writeUTF:
    file1.close()
if writeptsTF:
    file2.close()



# OLD DATA:
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








# rand(1)
# nPts = 108
# N = 51


# old: -6.0202e-19         new: 5.2352e-18         diff: 5.8372e-18        box: 2.0159e-09


# boxNew: 2.0250e-09       box: 2.0159e-09         dH: -6.5696e-21
# passed
# 0.7407407407407408 % passed     d 1.3005139389942589e-09
# 40.0 % passed   dv 5.90509347176418e-28

# N = 5001

# old: -1.0573e-18         new: -1.0553e-18        diff: 2.0661e-21        box: 2.2302e-09


# boxNew: 2.2295e-09       box: 2.2302e-09         dH:  3.3200e-21
# failed
# 49.351851851851855 % passed     d 1.4873996747615223e-10
# 50.0 % passed   dv 2.8064040336364697e-29

# Elapsed time: 73.0079 seconds
# 1.6425673736158061e-15 -1.0552607934327202e-18
