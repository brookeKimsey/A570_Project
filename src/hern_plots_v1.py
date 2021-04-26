### File to read in snapshots ###
### An attempt to write an hernquist 93
### replicate. Never used ###

import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import functions as fcn
import argparse

parser = argparse.ArgumentParser(description =
    'Read Gadget Snapshots & plot')

parser.add_argument('part',type = str, 
    help = 'Particles to plot; h=halo,d=disk,b=bulge')
parser.add_argument('path_to_files',type = str,  
    help = 'relative path to *.hdf5 files')
parser.add_argument('snap',nargs='+', type = str,
    help = 'snapshot numbers to plot (must be at least 3 numbers for each snap); --snap 000 100 987')


args = parser.parse_args()


file_list = []

for val in args.snap:
    fname = str(args.path_to_files) + 'snapshot_{}.hdf5'.format(val)
    file_list.append(fname)

if args.part == 'd':
    typ_choice = 'PartType2'
if args.part == 'h':
    typ_choice = 'PartType1'
if args.part == 'b':
    typ_choice = 'PartType3'

n_list = len(file_list)
if n_list in (1,2,3,4,5):
    n_rows,n_cols = 1,n_list
if n_list in (6,7,8,9,10):
    n_rows,n_cols = 2, int(ceil(n_list/2))
if n_list > 10:
    print('Too many Snapshots (10 or less)')
    print('Edit script if more are desired')
    exit()

fig,axes = plt.subplots(nrows=n_rows,ncols=n_cols)
i,k = 0,0

for filename in file_list:
    f = h5py.File(filename,'r')
    #file_list.append(f)

    for grp in list(f.keys()):
        if grp == type_choice:
            type_N = f['{}'.format(grp)]
            vel = type_N['Velocities']
            pos = type_N['Coordinates']
            ids = type_N['ParticleIDs']
            
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
