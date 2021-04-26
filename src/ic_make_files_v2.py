### Creates initial conditions for halo, disk, bulge ###
### uses exp_r_v2() the second version of the exp    ###
### disk and first version of bulge_r()              ###
### to be input for ./binary3                        ###

import numpy as np
import random
import functions as fcn

exec(open("values.py").read())

pts_disk = fcn.exp_r_v2(N=N_disk) ###only change is here compared to ic_make_files.py
pts_bulge = fcn.bulge_r(N=N_bulge,a = a_val/10)
pts_halo = fcn.halo_r(N=N_halo,a=a_val)
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
