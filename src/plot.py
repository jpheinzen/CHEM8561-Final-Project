import numpy as np
import matplotlib.pyplot as plt
from imports import *
import matplotlib.animation as animation
import sys


(U, resid) = importU()
(pts,stats) = importPts()
(arrLen,_,_) = np.shape(pts)

pts *= 1e9  # to get units of nm

"""Volume over time"""
# plt.plot(stats[:,0],stats[:,1]**(1/3))
# plt.show()


# plt.plot(pts[1,:,1],pts[1,:,2],'o')

"""Potential over time"""
# plt.semilogy(U[:,0],abs(U[:,1]),'o')
# plt.loglog(U[:,0],abs(U[:,1]),'o')
# plt.semilogy(U[:-1,0],resid,'o')
# plt.loglog(U[:-1,0],resid,'o')
# plt.show(block=False)
# plt.pause(10)

# sys.exit()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# Create a figure and axis objects
# fig, ax = plt.subplots(projection='3d')


i = 1
# ax.scatter(pts[i,:,0], pts[i,:,1], pts[i,:,2], marker='o')


ax.axis("tight")

# plt.show()

# The function to update the figure for each frame
def update(num):
    num *= 100
    global pts
    ax.clear()

    ax.scatter(pts[num,:,0], pts[num,:,1], pts[num,:,2], marker='o')

    plt.title(f'Configuration: {num+1}')
    mx = np.max(pts[num,:,:])
    ax.set(xlim=(0, mx), ylim=(0, mx), zlim=(0,mx))
    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=arrLen//100)

plt.show()


