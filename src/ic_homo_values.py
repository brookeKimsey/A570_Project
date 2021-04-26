### Create homogeneous disk, bulge, halo ###
### initial conditions and write to files ###
### to be input of ./binary3 

import numpy as np
import random
import functions as fcn

exec(open("values.py").read())

pts_disk = fcn.disk_homo_r(10,N=N_disk)
pts_bulge = fcn.sph_homo_r(4,N=N_bulge)
pts_halo = fcn.sph_homo_r(11,N=N_halo)
v_halo = fcn.vel(mu_halo,sigma_halo,N_halo)
v_disk = fcn.vel(mu_disk,sigma_disk,N_disk)
v_bulge = fcn.vel(mu_bulge,sigma_bulge,N_bulge) 

pts_all = np.vstack((pts_halo,pts_disk,pts_bulge))
pts_hd = np.vstack((pts_halo,pts_disk))
pts_db = np.vstack((pts_disk,pts_bulge))
np.savetxt('../attempt/d_xyz.txt',pts_disk)
np.savetxt('../attempt/all_xyz.txt',pts_all)
np.savetxt('../attempt/hd_xyz.txt',pts_hd)
np.savetxt('../attempt/db_xyz.txt',pts_db)

v_all = np.vstack((v_halo,v_disk,v_bulge))
v_hd = np.vstack((v_halo,v_disk))
v_db = np.vstack((v_disk,v_bulge))
np.savetxt('../attempt/d_vel.txt',v_disk)
np.savetxt('../attempt/all_vel.txt',v_all)
np.savetxt('../attempt/hd_vel.txt',v_hd)
np.savetxt('../attempt/db_vel.txt',v_db)

print('N Halo = {},N disk = {}, N Bulge = {}'.format(N_halo,N_disk,N_bulge))
