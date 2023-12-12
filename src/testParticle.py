import numpy as np
# import matplotlib.pyplot as plt
from imports import *
from functions import *
import sys

# nums = rng.random([3,nPts])





(U, resid) = importU()
(pts,stats) = importPts()
(arrLen,nPts,_) = np.shape(pts)


print(arrLen,nPts)

# Stats:
# iter, V, dp, dv



bg = np.random.MT19937(2)
rng = np.random.Generator(bg)

V = stats[:,1]
boxMat = V ** (1/3)


print("all points are in box?",checkInBox(pts,boxMat,arrLen))


epsNAP = 353.2
sigmaNAP = 4.7290e-10

num_TP_Pts = 5


i = 0

TP_pts = generateConfig(num_TP_Pts, boxMat[i], rng)
print(TP_pts)

j = 0
(LJ12,LJ6) = getPotential(TP_pts,j,boxMat[i])
U = addLJArr(LJ12=LJ12,LJ6=LJ6)

print(U)











