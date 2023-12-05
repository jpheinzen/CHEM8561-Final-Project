import numpy as np
import matplotlib.pyplot as plt
import sys
# from matplotlib.pyplot import plot, ion, show, draw

# interactive mode
# ion()

k =  1.380649e-23        # J/K.
sigma = 3.831e-10**2
eps = 204.68*k

def generateConfig(nPts: int, box: float, rng) -> np.ndarray:
    # rng = np.random.default_rng()  # the simplest way to generate random numbers
    nums = rng.random([3,nPts])
    nums = box*nums
    return nums

def plotParticles(pts: np.ndarray, xmax: float = 0, ymax: float = 0) -> None:
    plt.plot(pts[1,:],pts[2,:],'o')
    if xmax != 0:
        plt.xlim(0,xmax)
    if ymax != 0:
        plt.ylim(0,ymax)
    plt.show()

def calcLJPotential(sqDistance: np.ndarray) -> "tuple[np.ndarray,np.ndarray]":
    global k, sigma, eps 
    
    v = (sigma/sqDistance)**3    # type: ignore
    # LJ = 4*eps*(v**2 - v)           # type: ignore
    LJ12 = (v**2)             # type: ignore
    LJ6 = (v)                # type: ignore

    return LJ12,LJ6

def getPotential(pts: np.ndarray, i: int, box: float) -> "tuple[np.ndarray,np.ndarray]":
    hbox = box/2
    # U = 0.0
    # print(hbox)
    # print(pts)

    # print()
    # print(i)
    dr = np.concatenate((pts[:,:i],pts[:,i+1:]),1) - pts[:,i:i+1]
    # print(dr)
    # dr = np.abs(dr)

    dr[dr > hbox] -= box
    dr[dr < -hbox] += box
    # dr = np.where(dr > hbox, dr-box, dr)
    # dr = np.where(dr < -hbox, dr+box, dr)
    # print(dr)

    # sqDistance2 = dr[:,0]
    sqDistance = dr[0:1,:]**2 + dr[1:2,:]**2 +dr[2:3,:]**2

    # print(sqDistance2)
    # print(sqDistance)
    # print(sqDistance2+sqDistance)
    
    return calcLJPotential(sqDistance)

def getPotentials(pts: np.ndarray, nPts: int, box: float) -> "tuple[np.ndarray,np.ndarray]":
    hbox = box/2
    # U = 0.0
    # print(hbox)
    # print(pts)

    # Uarr = np.zeros(nPts)
    LJ12 = np.zeros(hashF(nPts))
    LJ6 = np.zeros(hashF(nPts))

    for i in range(1,nPts):
        # print()
        # print(i)
        hfi = hashF(i)
        hfip1 = hashF(i+1)
        # print(hfi,hfip1)
        dr = pts[:,:i] - pts[:,i:i+1]
        # print(dr)
        # dr = np.abs(dr)

        dr[dr > hbox] -= box
        dr[dr < -hbox] += box
        # dr = np.where(dr > hbox, dr-box, dr)
        # dr = np.where(dr < -hbox, dr+box, dr)
        # print(dr)

        # sqDistance2 = dr[:,0]
        sqDistance = dr[0:1,:]**2 + dr[1:2,:]**2 +dr[2:3,:]**2

        # print(sqDistance2)
        # print(sqDistance)
        # print(sqDistance2+sqDistance)
        
        (LJ12[hfi:hfip1],LJ6[hfi:hfip1]) = calcLJPotential(sqDistance)
        # print(pots)
        # Uarr[i] = np.sum(LJ12[hfi:hfip1]+LJ6[hfi:hfip1])
        # print(U)
        # print(i)y
    
    # U = np.sum(Uarr)

    return LJ12,LJ6

def calcLJPotentialOld(sqDistance: np.ndarray) -> "tuple[np.ndarray,np.ndarray]":
    global k, sigma, eps 
    
    v = (sigma/sqDistance)**3    # type: ignore
    # LJ = 4*eps*(v**2 - v)           # type: ignore
    LJ12 = 4*eps*(v**2)             # type: ignore
    LJ6 = -4*eps*(v)                # type: ignore

    return LJ12,LJ6

def getPotentialOld(pts: np.ndarray, nPts: int, box: float) -> float:
    hbox = box/2
    U = 0.0
    # print(hbox)
    # print(pts)

    for i in range(nPts-1):
        # print()
        # print(i)
        dr = pts[:,i+1:] - pts[:,i:i+1]
        # print(dr)
        # dr = np.abs(dr)

        dr[dr > hbox] -= box
        dr[dr < -hbox] += box
        # dr = np.where(dr > hbox, dr-box, dr)
        # dr = np.where(dr < -hbox, dr+box, dr)
        # print(dr)

        # sqDistance2 = dr[:,0]
        sqDistance = dr[0:1,:]**2 + dr[1:2,:]**2 +dr[2:3,:]**2

        # print(sqDistance2)
        # print(sqDistance)
        # print(sqDistance2+sqDistance)
        
        (LJ12,LJ6) = calcLJPotentialOld(sqDistance)
        # print(pots)
        U += np.sum(LJ12+LJ6)
        # print(U)
        # print(i)

    return U

def acceptMove(dPotential: float, beta: float, rng) -> bool:

    p = np.min((1,np.exp(-beta*dPotential))) # type: ignore
    r = rng.random()
    # print('random number',r,'\tvalue',np.exp(-beta*dPotential))
    return r < p

def hashF(i: int) -> int:
    return i*(i-1)//2

def integers(a, b):
    return np.arange(a,b)

def getInds(i: int, N: int) -> np.ndarray:
    return np.concatenate((integers(hashF(i),hashF(i+1)), hashF(integers(i+1,N))+i)) # type: ignore

def check(a,b,tag,exit=True, tol = 1e-15):
    if (abs((a-b)/a) > tol).any() and (abs(a-b) > 1e-32).all():
        print(tag,'FAILED','a',a,'b',b,(a-b)/a,abs(a-b))
        if exit:
            sys.exit()
        else:
            print()
            # time.sleep(0.1)

if __name__ == "__main__":
    # rng = np.random.default_rng()
    bg = np.random.MT19937(1)
    rng = np.random.Generator(bg)
    pts = generateConfig(4,10, rng)
    # print(nums)

    # print(np.min(nums))
    # print(np.max(nums))
    # print(nums)
    # print(nums[1:,:])
    # print(np.max((1,np.exp(1))))
    
    # print()




    # a = np.arange(6)
    # a1 = a**2
    # a2 = a**3

    # a +=1
    # print(a,a1,a2)

    # print(a1[getInds(2,4)])

    # b = np.zeros(4)
    # print(b)

    # A = np.arange(7)
    # print(A)
    # A = A**2
    # print(A)

    

    # minval = 2**(1/6)*sigma

    # r = np.linspace(minval/1.1,3*minval,1000)**2    # type: ignore
    # LJ = calcLJPotential(r)

    # print(LJ)

    # plt.plot(r/sigma,LJ/eps)
    # plt.plot(r,LJ)
    # plt.ylim(-2,4)
    # draw()
    # plt.plot(r,LJ)
    # draw()
    # plt.show(block=False)

    # plt.plot(r,LJ)
    # plt.show(block=False)


    # print(minval)
    # plt.show()
    

    # bg = np.random.MT19937(0)
    # rg = np.random.Generator(bg)
    # print(rg.random())




    i = 2
    print(pts[:,:i])
    print(pts[:,i+1:])
    a = np.concatenate((pts[:,:i],pts[:,i+1:]),1)
    print(a)















