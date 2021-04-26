### Functions        ###
### Original version ###

import numpy as np
import random

#M_disk = 5
#M_halo = 12*M_disk
#M_bulge = M_disk/10
#N_disk = 5000
#N_halo = int(3*N_disk)
#N_bulge = int(N_disk/10)
#a_val = 1.0
#h_val = 1.0
#z_0 = 1.0
#pb0,pd0,ph0 = 0.427,100.,0.711
#mu_disk, sigma_disk = 1.0, 0.1
#mu_halo, sigma_halo = 1.0, 0.1
#mu_bulge, sigma_bulge = 1.0, 0.1

exec(open("values.py").read())

### Velocities ###

def vel(mu,sig,N):
    vx = np.array([random.gauss(mu,sig) for _ in range(N)])
    vy = np.array([random.gauss(mu,sig) for _ in range(N)])
    vz = np.array([random.gauss(mu,sig) for _ in range(N)])
    v = np.stack((vx,vy,vz),axis=-1)
    return v

### Positions ###

def sph_homo_r(R,N):
    '''An homogeneous sphere '''
    u = np.array([random.uniform(0,1) for _ in range(N)])
    v = np.array([random.uniform(0,1) for _ in range(N)])
    w = R*np.array([random.uniform(0,1) for _ in range(N)])**(1/3)

    theta = 2*np.pi*np.array(u)
    phi = np.arccos(2*np.array(v) -1)

    x = w*np.sin(phi)*np.cos(theta)
    y = w*np.sin(phi)*np.sin(theta)
    z = w*np.cos(phi)

    p = np.stack((x,y,z),axis=-1)

    return p  
 
def disk_homo_r(R,N):
    '''An homogeneous disk '''
    u = np.array([random.uniform(0,1) for _ in range(N)])
    w = R*np.array([random.uniform(0,1) for _ in range(N)])**(1/2)

    theta = 2*np.pi*np.array(u)

    x = w*np.cos(theta)
    y = w*np.sin(theta)
    z = np.array([random.uniform(-1,1) for _ in range(N)])

    p = np.stack((x,y,z),axis=-1)

    return p

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

### Density Checks ###
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

def vel_disp(X,vel,N):
    ndim = X.shape[1]
    if ndim == 3:
        r = np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    if ndim == 2:
        r = np.sqrt(X[:,0]**2 + X[:,1]**2)

    #print(r)
    #print(vel)
    comb = np.stack((r,vel),axis=-1)
    #print(comb)
    comb_sort = comb[np.argsort(comb[:,0])]
    #print(comb_sort)

    N_velArr = np.array([])
    r_MedArr = np.array([])

    ### Ethan's Help ###
    for i in range(len(comb_sort[::N])):
        i_inr = i*(N-1)
        i_outr = (i+1)*(N-1)
        #print(i_inr,i_outr)
        r_inr =  comb_sort[i_inr,0]
        r_outr = comb_sort[i_outr,0]
        #print(r_inr,r_outr)
        vel_in_shell = comb_sort[i_inr:(i_outr+1),1]
        #print(vel_in_shell)
        #v_inr = v_sort[i*(N-1)]
        #v_outr = v_sort[(i+1)*(N-1)]

        v_disp = np.std(vel_in_shell)
        #print(v_disp)
        #if ndim == 3:
        #V = (4/3)*np.pi*(r_outr**3 - r_inr**3)

        #N_den = N/V
        r_med = np.median(comb_sort[i_inr:(i_outr+1),0])
        #print(r_med)

        N_velArr = np.append(N_velArr, v_disp)
        r_MedArr = np.append(r_MedArr, r_med)

    return N_velArr, r_MedArr
