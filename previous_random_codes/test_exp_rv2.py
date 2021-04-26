###test 2nd version of exponential disk ###

import numpy as np
import functions as fcn
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pts = fcn.exp_r_v2(1000,10)

fig,ax = plt.subplots()

ax.scatter(pts[:,0],pts[:,1])

fig3d = plt.figure()
ax3d = fig3d.gca(projection='3d')
ax3d.scatter(pts[:,0],pts[:,1],pts[:,2])
ax3d.set(xlim=(-30,30),ylim=(-30,30),zlim=(-30,30))

plt.show()
