### File to read in snapshots ###
### Earlier version

import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt

path_to_h5_files = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/attempt/all_results'
all_files = glob.glob(path_to_h5_files + "/*.hdf5")

file_list = []

for filename in all_files:
    f = h5py.File(filename,'r')
    file_list.append(f)
    #print(filename)
    #print(list(f.keys()))
    for grp in list(f.keys()):
        if grp != 'Header':
            type_N = f['{}'.format(grp)]
            vel = type_N['Velocities']
            pos = type_N['Coordinates']
            ids = type_N['ParticleIDs']
            #print(grp,type_N,vel[:,0])
    #header = f['Header']
    #type2 = f['PartType2']
    #vel = type2['Velocities']
    #pos = type2['Coordinates']
    #ids = type2['ParticleIDs']
    #print(header.items())
    #print(type2.keys())
    #print(ids)

            fig,axes = plt.subplots(figsize=(20,20),nrows=2,ncols=2,sharey=True,sharex=True)
            plt.subplots_adjust(hspace=0,wspace=0)
            fig.suptitle('Snapshot={},{} Only'.format(filename,type_N))
            axes[0,0].scatter(pos[:,0],pos[:,2],s=1,c='k')
            axes[1,0].scatter(pos[:,0],pos[:,1],s=1,c='k')
            axes[0,1].scatter(pos[:,1],pos[:,2],s=1,c='k')
            axes[1,0].set_xlabel('x',fontsize=12)
            axes[1,0].set_ylabel('y')
            axes[0,0].set_ylabel('z')
            axes[0,1].set_xlabel('y')
            axes[1,1].axis('off')

            plt.savefig('{}_{}_pos.png'.format(filename,grp))
            #plt.show()
            plt.close(fig)
