### testing creating input for ./binary

import numpy as np
import random

M_disk = 5
M_halo = 12*M_disk
M_bulge = M_disk/10
N_disk = 10000
N_halo = int(3*N_disk)
N_bulge = int(N_disk/10)
a_val = 1.0
h_val = 1.0
z_0 = 1.0
pb0,pd0,ph0 = 0.427,100.,0.711
mu_disk, sigma_disk = 1.0, 0.1
mu_halo, sigma_halo = 1.0, 0.1
mu_bulge, sigma_bulge = 1.0, 0.1

def vel(mu,sig,N):
    vx = np.array([random.gauss(mu,sig) for _ in range(N)])
    vy = np.array([random.gauss(mu,sig) for _ in range(N)])
    vz = np.array([random.gauss(mu,sig) for _ in range(N)])
    v = np.stack((vx,vy,vz),axis=-1)
    return v

def halo_r(N=10,a=1.0):
    '''See write up for derivation'''
    u = np.array([random.uniform(0,1) for _ in range(N)])
    v = np.array([random.uniform(0,1) for _ in range(N)])
    q_r = np.array([random.uniform(0,1) for _ in range(N)])

    theta = 2*np.pi*u
    phi = np.arccos(2*v - 1)

    w = a*(np.sqrt(q_r)/(1-np.sqrt(q_r)))

    x = w*np.sin(phi)*np.cos(theta)
    y = w*np.sin(phi)*np.sin(theta)
    z = w*np.cos(phi)

    p = np.stack((x,y,z),axis=-1)

    return p

def exp_r(N,rd=2.0):
    q_r = np.random.uniform(0,1,N)
    u = np.random.uniform(0,1,N)
    R = np.linspace(0,10*rd,N)

    theta = 2*np.pi*u
    LHS = (np.exp(-R/rd)*(R/rd + 1))
    RHS = 1-q_r
    t = np.stack((LHS,RHS),axis=-1)
    t2 = t[t[:,0].argsort()]
    w = np.interp(R, t2[:,0], t2[:,1])
    x = w*np.cos(theta)
    y = w*np.sin(theta)
    z = np.random.normal(0,0.1,N)

    p3d = np.stack((x,y,z),axis=-1)

    return p3d

def bulge_r(N=10,a=1.0,p0=1.0,Mtot=1.0):
    '''See write up for derivation'''
    q_r = np.random.uniform(0,1,N)
    u = np.random.uniform(0,1,N)
    v = np.random.uniform(0,1,N)

    theta = 2*np.pi*u
    phi = np.arccos(2*v - 1)

    c = q_r*2*np.pi*a**3*p0/Mtot
    w = (np.sqrt(c) - c)/(c-1)


    x = w*np.sin(phi)*np.cos(theta)
    y = w*np.sin(phi)*np.sin(theta)
    z = w*np.cos(phi)

    p = np.stack((x,y,z),axis=-1)

    return p

pts_disk = exp_r(N=N_disk)
pts_bulge = halo_r(N=N_bulge,a = a_val/10)
pts_halo = halo_r(N=N_halo,a=a_val)
v_halo = vel(mu_halo,sigma_halo,N_halo)
v_disk = vel(mu_disk,sigma_disk,N_disk)
v_bulge = vel(mu_bulge,sigma_bulge,N_bulge) 

np.savetxt('disk_xyz.txt',pts_disk)
np.savetxt('bulge_xyz.txt',pts_bulge)
np.savetxt('halo_xyz.txt',pts_halo)
np.savetxt('disk_vel.txt',v_disk)
np.savetxt('bulge_vel.txt',v_bulge)
np.savetxt('halo_vel.txt',v_halo)

print('N Disk = {},N Halo = {}, N Bulge = {}'.format(N_disk,N_halo,N_bulge))
