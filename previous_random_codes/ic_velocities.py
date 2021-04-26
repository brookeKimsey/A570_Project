### original attempt at creating ic_velocities
### before creating a function to do this
### create hdf5 file of ic velocities
### bulge, halo, disk

import numpy as np
import h5py

N_disk, N_halo, N_bulge = 10,10,10
mu_disk, sigma_disk = 1.0, 0.1
mu_halo, sigma_halo = 1.0, 0.1
mu_bulge, sigma_bulge = 1.0, 0.1

vx_disk = np.random.normal(mu_disk,sigma_disk,N_disk)
vy_disk = np.random.normal(mu_disk,sigma_disk,N_disk)
vz_disk = np.random.normal(mu_disk,sigma_disk,N_disk)

v_disk = np.stack((vx_disk,vy_disk,vz_disk),axis=-1)

vx_halo = np.random.normal(mu_halo,sigma_halo,N_halo)
vy_halo = np.random.normal(mu_halo,sigma_halo,N_halo)
vz_halo = np.random.normal(mu_halo,sigma_halo,N_halo)

v_halo = np.stack((vx_halo,vy_halo,vz_halo),axis=-1)

vx_bulge = np.random.normal(mu_bulge,sigma_bulge,N_bulge)
vy_bulge = np.random.normal(mu_bulge,sigma_bulge,N_bulge)
vz_bulge = np.random.normal(mu_bulge,sigma_bulge,N_bulge)

v_bulge = np.stack((vx_bulge,vy_bulge,vz_bulge),axis=-1)

print(v_disk,v_halo,v_bulge)

