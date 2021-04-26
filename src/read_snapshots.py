### File to read in snapshots          ###
### Takes Gadget2 output (.hdf5 files) ###
### and creates plots modeling all     ###
### particles in a given simulation,   ###
### along with density distribution and###
### velocity dispersion                ###

import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import functions as fcn
import argparse

parser = argparse.ArgumentParser(description =
    'Read Gadget Snapshots & plot')

parser.add_argument('combo',type = str, 
    help = 'Halo,disk, bulge = all; halo,disk=hd; disk,bulge = bd; disk= d')
parser.add_argument('--N',type = int, default=10, # in homoI and homoII, the value is 10
    help = 'Number in each volume slice; default=10')

args = parser.parse_args()
N_val = args.N

path_to_h5_files = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/attempt/{}_results'.format(args.combo)
all_files = glob.glob(path_to_h5_files + "/*.hdf5")

file_list = []

for filename in all_files:
    f = h5py.File(filename,'r')
    file_list.append(f)

    v_sys = np.array([0,0,0])
    p_sys = np.array([0,0,0])

    fig,axes = plt.subplots(figsize=(14,8),nrows=2,ncols=3)
    for grp in list(f.keys()):
        if grp != 'Header':
            type_N = f['{}'.format(grp)]
            vel = type_N['Velocities']
            pos = type_N['Coordinates']
            ids = type_N['ParticleIDs']
            
            ### System ### 
            v_sys = np.vstack((v_sys,vel))
            p_sys = np.vstack((p_sys,pos))
            
            ### Density Distribution ###
            dens,r_dens = fcn.density(pos,N_val)
            v_disp, r_disp = fcn.vel_disp(pos,vel,N_val)

            if grp == 'PartType1':
                c = 'k'
                zo = 100
            elif grp == 'PartType2':
                c = 'r'
                zo = 1000
            elif grp == 'PartType3':
                c = 'b'
                zo = 100000 ###this was one in homoI & homoII
             
            #star positions
            axes[0,0].scatter(pos[:,0],pos[:,2],s=1,c=c,zorder=zo)
            axes[1,0].scatter(pos[:,0],pos[:,1],s=1,c=c,zorder=zo)
            axes[0,1].scatter(pos[:,1],pos[:,2],s=1,c=c,zorder=zo)
            #density distribution
            axes[1,1].plot(r_dens,dens,color=c,zorder=zo)
            #velocity dispersion
            axes[0,2].plot(r_disp,v_disp,color=c,zorder=zo)            
   
    ### System ###
    vsys_disp, rsys_disp = fcn.vel_disp(p_sys[1:],v_sys[1:],N_val)
    denssys,rsys_dens = fcn.density(p_sys[1:],N_val)
    axes[1,1].plot(rsys_dens,denssys,color='magenta',zorder=0)
    axes[0,2].plot(rsys_disp,vsys_disp,color='magenta',zorder=0)

    axes[0,1].set_title('Snapshot={}'.format(filename[-10:-4]))
    axes[1,0].set_xlabel('x')
    axes[1,0].set_ylabel('y')
    axes[0,0].set_ylabel('z')
    axes[0,1].set_xlabel('y')#,fontsize=20)
    axes[1,1].set_ylabel('Density Distribution')#,fontsize=20)
    axes[1,1].set_xlabel('Radius')#,fontsize=20)
    axes[0,2].set_ylabel(r'v$_r$ dispersion')
    axes[0,2].set_xlabel('Radius')
    axes[1,2].axis('off')
    axes[1,2].scatter([],[],c='k',label='Halo')
    axes[1,2].scatter([],[],c='r',label='Disk')
    axes[1,2].scatter([],[],c='b',label='Bulge') 
    axes[1,2].scatter([],[],c='magenta',label='System')

    axes[0,0].set(xlim=(-10,10),ylim=(-10,10))
    axes[1,0].set(xlim=(-10,10),ylim=(-10,10))
    axes[0,1].set(xlim=(-10,10),ylim=(-10,10))
    axes[1,1].set(xlim=(0,10),ylim=(0,1000))### in certain models, the limit is 100
    axes[0,2].set(xlim=(0,10),ylim=(0,100)) ### in certain models, the limit is 100
    
    axes[1,2].legend()
    plt.tight_layout()
    plt.savefig('../attempt/{}_results/plots/snap{}_pos.png'.format(args.combo,filename[-10:-5]))
    plt.close(fig)
