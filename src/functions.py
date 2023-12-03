import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.pyplot import plot, ion, show, draw

# interactive mode
# ion()

k =  1.380649e-23        # J/K.
sigma = 3.831e-10
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

def calcLJPotential(sqDistance: np.ndarray) -> np.ndarray:
    global k, sigma, eps 
    
    v = (sigma**2/sqDistance)**3    # type: ignore
    LJ = 4*eps*(v**2 - v)           # type: ignore
    # LJ6 = 4*eps*(v**12 - v**6) # type: ignore

    return LJ

def getPotential(pts: np.ndarray, nPts: int, box: float) -> float:
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
        
        pots = calcLJPotential(sqDistance)
        # print(pots)
        U += np.sum(pots)
        # print(U)
    
    return U

def acceptMove(dPotential: float, beta: float, rng) -> bool:

    p = np.min((1,np.exp(-beta*dPotential))) # type: ignore
    r = rng.random()
    # print('random number',r,'\tvalue',np.exp(-beta*dPotential))
    return r < p


if __name__ == "__main__":
    rng = np.random.default_rng()
    nums = generateConfig(4,100, rng)
    # print(nums)

    # print(np.min(nums))
    # print(np.max(nums))
    # print(nums)
    # print(nums[1:,:])
    # print(np.max((1,np.exp(1))))
    
    # print()

    b = np.zeros(4)
    print(b)

    A = np.arange(7)
    print(A)
    A = A**2
    print(A)

    

    minval = 2**(1/6)*sigma

    r = np.linspace(minval/1.1,3*minval,1000)**2    # type: ignore
    LJ = calcLJPotential(r)

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


    print(minval)
    # plt.show()
    

    bg = np.random.MT19937(0)
    rg = np.random.Generator(bg)
    print(rg.random())

















