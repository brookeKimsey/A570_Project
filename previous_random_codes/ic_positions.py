### Another attempt at creating Gadget2 input
### create hdf5 file of ic positions
### bulge, halo, disk

import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

M_disk = 5
M_halo = 12*M_disk
M_bulge = M_disk/10
N_disk = 100000
N_halo = int(3*N_disk)
N_bulge = int(N_disk/10)
a_val = 1.0
h_val = 1.0
z_0 = 1.0
pb0,pd0,ph0 = 0.427,100.,0.711

def halo_r(N=10,a=1.0):
    '''See write up for derivation'''
    u = np.array([random.uniform(0,1) for _ in range(N)])
    v = np.array([random.uniform(0,1) for _ in range(N)])
    q_r = np.array([random.uniform(0,1) for _ in range(N)])
    #q_r = np.random.uniform(0,1,N)
    #u = np.random.uniform(0,1,N)
    #v = np.random.uniform(0,1,N)
    
    theta = 2*np.pi*u
    phi = np.arccos(2*v - 1)

    w = a*(np.sqrt(q_r)/(1-np.sqrt(q_r)))
    
    x = w*np.sin(phi)*np.cos(theta)
    y = w*np.sin(phi)*np.sin(theta)
    z = w*np.cos(phi)

    p = np.stack((x,y,z),axis=-1)

    return p

def disk_r(N=10,h=1.0,z0=1.0,p0=1.0,Mtot=1.0):
    """ See write up for derivation """
    r = np.linspace(0,100*h,N)
    z = np.array([random.uniform(0,z0) for _ in range(N)])
    #z = np.random.uniform(0,100*z0,N)
    #q_r = np.random.uniform(0,1,N)
    q_r = np.array([random.uniform(0,1) for _ in range(N)]) 

    LHS = np.exp(-r/h)*(1+ r/h)
    RHS = 1 - ( (2*np.pi*q_r*p0*h**2)/(Mtot*np.tanh(z/z0)) )
   
    w = np.interp(r,LHS,RHS)
    u = np.array([random.uniform(0,1) for _ in range(N)])
    #u = np.random.uniform(0,1,N)
    theta = 2*np.pi*u

    x = w*np.cos(theta)
    y = w*np.sin(theta)

    p = np.stack((x,y,z),axis=-1)

    return p
   
def exp_r(N,rd=2.0):
    q_r = np.random.uniform(0,1,N)
    u = np.random.uniform(0,1,N)
    R = np.linspace(0,10*rd,N)
     
    theta = 2*np.pi*u
    LHS = (np.exp(-R/rd)*(R/rd + 1))
    RHS = 1-q_r
    #print(LHS)
    t = np.stack((LHS,RHS),axis=-1)
    t2 = t[t[:,0].argsort()]
    w = np.interp(R, t2[:,0], t2[:,1])
    #w = np.interp(R,LHS,RHS)
    x = w*np.cos(theta)
    y = w*np.sin(theta)
    z = np.random.normal(0,0.1,N)
 
    p2d = np.stack((x,y),axis=-1)
    p3d = np.stack((x,y,z),axis=-1)

    return p2d,p3d
         
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

def halo_rho(r,a):
    '''Hernquist Model BT2008 Equation 2.64'''
    rho_h = ph0/( (r/a)**1 * (1 + r/a)**3 )
    return rho_h    

def exp_rho(r,rd=2.0,p0=10.0):
    ''' 2d exponential disk '''
    rho_2d = p0*np.exp(-r/rd)
    return rho_2d


def disk_rho(r_flat,z,h=1.0,z0=1.0,Mtot=1.0):
    ''' Equation 2.1 Hernquist 1993'''
    x,y,z = r[:,0],r[:,1],r[:,2]
    r_flat = ( x**2 + y**2)**0.5
    rho_d = (Mtot/ (4*np.pi*h**2*z0) ) * np.exp(-r_flat/h)*(np.cosh(z/z0))**(-2)
    return rho_d 

def bulge_rho(r,a=1.0,c=1.0,Mtot=1.0):
    ''' Equation 2.6, 2.7 Hernquist 1993,assumes a=c'''
    #x, y, z = r[:,0],r[:,1], r[:,2]
    #m_val = (((x**2 + y**2) / a**2) + (z/c)**2)**0.5
    m_val = r/a
    rho_b = (Mtot/(2*np.pi*a*c**2) )*(1/(m_val * (1+m_val)**3))
    return rho_b

def plummer_rho(r,a=1.0,Mtot=1.0):
    rho = ( 3*Mtot/(4*np.pi*a**3) )*(1+(r/a)**2)**(-5/2)
    return rho

def density(X,N):
    ndim = X.shape[1]
    if ndim == 3:
        r = np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    if ndim == 2:
        r = np.sqrt(X[:,0]**2 + X[:,1]**2)
    r_sort = np.sort(r)

    N_denArr = np.array([])
    r_MedArr = np.array([])

    ### Ethan's Help ###
    for i in range(len(r_sort[::N])):
        r_inr =  r_sort[i*(N-1)]
        r_outr = r_sort[(i+1)*(N-1)]

        if ndim == 2:
            V = np.pi*(r_outr**2 - r_inr**2)
        if ndim == 3:
            V = (4/3)*np.pi*(r_outr**3 - r_inr**3)

        N_den = N/V
        r_med = np.median(r_sort[(i*(N-1)):((i+1) * (N-1))])

        N_denArr = np.append(N_denArr, N_den)
        r_MedArr = np.append(r_MedArr, r_med)

    return N_denArr, r_MedArr

### Test ###
pts_halo = halo_r(N=N_halo,a=a_val)
dens_halo, r_halo = density(pts_halo,100)

#r_xyz_halo = (np.sum((pts_halo)**2,axis=1))**0.5
#keep = np.where(r_xyz_halo <= 10.)
pts_disk2d,pts_disk3d = exp_r(N=N_disk) 
#pts_disk = disk_r(N=N_disk,h=h_val,z0=z_0,p0=pd0,Mtot=M_disk)
pts_bulge = halo_r(N=N_bulge,a=a_val/10)#,p0=pb0,Mtot=M_bulge)

#dens_disk, r_disk = density(pts_disk[:,0:2],10)
dens_disk, r_disk = density(pts_disk2d,100)
dens_bulge, r_bulge = density(pts_bulge,10)

fig3d = plt.figure()
ax3d = fig3d.add_subplot(111,projection = '3d')
ax3d.scatter(pts_disk3d[:,0],pts_disk3d[:,1],pts_disk3d[:,2],c='k',label='Disk')
ax3d.scatter(pts_halo[:,0],pts_halo[:,1],pts_halo[:,2],c='r',label='Halo')
ax3d.scatter(pts_bulge[:,0],pts_bulge[:,1],pts_bulge[:,2],c='g',label='Bulge')
#ax3d.set_xlim(-1,1)
#ax3d.set_ylim(-1,1)
#ax3d.set_zlim(-1,1)

#fig,(ax,ax1,ax2) = plt.subplots(nrows=1,ncols=3)#,sharey=True)
fig, ax = plt.subplots()
ax.set_xlabel(r'r = $\sqrt{x^2 + y^2 + z^2}$')
ax.set_ylabel(r'$\rho$(R)')

ax.scatter(r_halo, dens_halo, c='r',label='Simulated Halo')
#ax.plot(r_halo,halo_rho(r_halo,a_val), c= 'r',label='Hernquist Halo')

ax.scatter(r_disk,dens_disk, c='k',label='Simulated Disk')
#ax.plot(r_disk,exp_rho(r_disk),c='r',label='Exponential Disk 2d')
#ax1.plot(pts_disk,disk_rho(pts_disk,h_val,z_0),c = 'r',label='Hernquist Disk')

ax.scatter(r_bulge,dens_bulge, c='g',label='Simulated Bulge')
#ax.plot(r_bulge,bulge_rho(r_bulge,a_val,a_val,M_bulge),c='r',label='Hernquist Bulge')
#ax.plot(r_bulge,plummer_rho(r_bulge,a_val,M_bulge),c='g',label='Plummer Model')
ax.set_yscale('log')
ax.set_xscale('log')


ax.legend()
plt.show()

