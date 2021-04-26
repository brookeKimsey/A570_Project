### File to read in snapshots ###
### Second version of read in snapshots ###

import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import functions as fcn

path_to_h5_files = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/attempt/all_results'
all_files = glob.glob(path_to_h5_files + "/*.hdf5")

file_list = []

for filename in all_files:
    f = h5py.File(filename,'r')
    file_list.append(f)
    fig,axes = plt.subplots(figsize=(20,20),nrows=2,ncols=2,sharey=True,sharex=True)
    #print(filename)
    #print(list(f.keys()))
    for grp in list(f.keys()):
        print(grp)
        #fig,axes = plt.subplots(figsize=(20,20),nrows=2,ncols=2,sharey=True,sharex=True)
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

            if grp == 'PartType1':
                c = 'k'
            elif grp == 'PartType2':
                c = 'r'
            elif grp == 'PartType3':
                c = 'b'
             
            #fig,axes = plt.subplots(figsize=(20,20),nrows=2,ncols=2,sharey=True,sharex=True)
            plt.subplots_adjust(hspace=0,wspace=0)
            fig.suptitle('Snapshot={}'.format(filename[-10:-4]),fontsize=30)
            axes[0,0].scatter(pos[:,0],pos[:,2],s=5,c=c)
            axes[1,0].scatter(pos[:,0],pos[:,1],s=5,c=c)
            axes[0,1].scatter(pos[:,1],pos[:,2],s=5,c=c)
            axes[1,0].set_xlabel('x',fontsize=20)
            axes[1,0].set_ylabel('y',fontsize=20)
            axes[0,0].set_ylabel('z',fontsize=20)
            axes[0,1].set_xlabel('y',fontsize=20)
            axes[1,1].axis('off')
            axes[0,0].set(xlim=(-10,10),ylim=(-10,10))
            axes[1,0].set(xlim=(-10,10),ylim=(-10,10))
            axes[0,1].set(xlim=(-10,10),ylim=(-10,10))

    plt.savefig('../attempt/all_results/plot_test/{}_pos.png'.format(filename))
            #plt.show()
    plt.close(fig)
