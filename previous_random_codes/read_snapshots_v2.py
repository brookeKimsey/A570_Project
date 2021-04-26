### File to read in snapshots ###
### Another older version of read in snapshots ###

import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import functions as fcn

path_to_h5_files = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/attempt/db_results'
all_files = glob.glob(path_to_h5_files + "/*.hdf5")

file_list = []

for filename in all_files:
    f = h5py.File(filename,'r')
    file_list.append(f)
    fig,axes = plt.subplots(figsize=(14,8),nrows=2,ncols=3)#,sharey=True,sharex=True)
    for grp in list(f.keys()):
        if grp != 'Header':
            type_N = f['{}'.format(grp)]
            vel = type_N['Velocities']
            pos = type_N['Coordinates']
            ids = type_N['ParticleIDs']

            ### Density Distribution ###
            dens,r_dens = fcn.density(pos,10)

            if grp == 'PartType1':
                c = 'k'
            elif grp == 'PartType2':
                c = 'r'
            elif grp == 'PartType3':
                c = 'b'
             
            #plt.subplots_adjust(hspace=0,wspace=0)
            axes[0,1].set_title('Snapshot={}'.format(filename[-10:-4]))#,fontsize=30)
            axes[0,0].scatter(pos[:,0],pos[:,2],s=1,c=c)
            axes[1,0].scatter(pos[:,0],pos[:,1],s=1,c=c)
            axes[0,1].scatter(pos[:,1],pos[:,2],s=1,c=c)
            axes[1,0].set_xlabel('x')#,fontsize=20)
            axes[1,0].set_ylabel('y')#,fontsize=20)
            axes[0,0].set_ylabel('z')#,fontsize=20)
            axes[0,1].set_xlabel('y')#,fontsize=20)
            axes[1,1].axis('off')
            axes[1,2].axis('off')
            axes[0,2].set_ylabel('Density Distribution')#,fontsize=20)
            axes[0,2].set_xlabel('Radius [Natural Units]')#,fontsize=20)
            axes[0,2].plot(r_dens,dens,lw = 5,color=c)
            axes[0,0].set(xlim=(-10,10),ylim=(-10,10))
            axes[1,0].set(xlim=(-10,10),ylim=(-10,10))
            axes[0,1].set(xlim=(-10,10),ylim=(-10,10))

    plt.tight_layout()
    plt.savefig('../attempt/db_results/plots/snap{}_pos.png'.format(filename[-10:-5]))
    plt.close(fig)
