import numpy as np
# import matplotlib.pyplot as plt
from imports import *
from functions import *
import sys

# nums = rng.random([3,nPts])


"""import U and Pts arrays from file. Will have to specify file path"""
(U, resid) = importU()
(pts,stats) = importPts()
# (U, resid) = importU(fileName='Us.txt')
# (pts,stats) = importPts(fileName='ptss.txt')
(arrLen,nPts,_) = np.shape(pts)
# Stats:
# iter, V, dp, dv

print(arrLen,nPts)


# rng = np.random.default_rng()
bg = np.random.MT19937(2)
rng = np.random.Generator(bg)

VMat = stats[:,1]
boxMat = VMat ** (1/3)


print("all points are in box?",checkInBox(pts,boxMat,arrLen))

T = 320
beta = 1/(k*T)

epsNAP = 353.2
sigmaNAP = 4.7290e-10

num_TP_Pts = 200
# num_TP_Pts = 500

num_configs = 1000

# These are from Eq 5 in the paper
numerator = np.zeros((num_configs,num_TP_Pts))
denominator = np.zeros((num_configs,num_TP_Pts))


# low = 5000//100-1
low = arrLen-num_configs-1
# low = 10
# for iConfig in range()
i = 0
for iConfig in range(low,low+num_configs):
    TP_pts = generateConfig(num_TP_Pts, boxMat[iConfig], rng)
    # print(TP_pts)

    psi = np.zeros(num_TP_Pts)
    for jTP in range(num_TP_Pts):
        # print(np.shape(TP_pts))
        # print(np.shape(np.transpose(pts[iConfig,:,:])))
        # print(np.shape(TP_pts[:,jTP:jTP+1]))

        (LJ12,LJ6) = getPotential(np.concatenate((TP_pts[:,jTP:jTP+1],np.transpose(pts[iConfig,:,:])),1),0,boxMat[iConfig],sigma=sigmaNAP)
        psi[jTP] = addLJArr(LJ12=LJ12,LJ6=LJ6,eps=epsNAP)
    # print('1',np.shape(psi))
    # print('2',np.shape(numerator[i,:]))

    if i % (num_configs//10) == 0:
        print('done with',i,iConfig)
    
    numerator[i,:] = VMat[iConfig] * np.exp(-beta*psi)
    denominator[i,:] = VMat[iConfig]
    i+=1
# print(psi)

# print(beta)
# print(-beta*psi)
# print(np.exp(-beta*psi))
# print(VMat[iConfig]*(1e9)**3)
# print((VMat[iConfig]*(1e9)**3)*np.exp(-beta*psi))

print('done')

num = np.average(numerator)
den = np.average(denominator)

val = -np.log(num/den)

print(num,den,val)

print(VMat)





