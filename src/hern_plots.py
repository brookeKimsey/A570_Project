### Recreate Hernquist 93 plots ###
### Compare isolated disk & disk with all 3
### components
### Pick 8 snapshots to compare times to ###
### pick which component to compare (disk,halo, bulge)

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse

snaps = ['000','032','053','100','200','300','500','1000']

fig, axes = plt.subplots(figsize=(10,10),nrows = 4, ncols=4,subplot_kw={'xticks':[],'yticks':[]})

attemptN = 16
i,j = 0,0
k = i+2

for val in snaps:
    fname_all = '../previous_attempts/attempt{}/all_results/snapshot_{}.hdf5'.format(attemptN,val)
    fname_iso = '../previous_attempts/attempt{}/d_results/snapshot_{}.hdf5'.format(attemptN,val)
    fall = h5py.File(fname_all,'r')
    fiso = h5py.File(fname_iso,'r')

    if len(val) == 3:
        t = '{}.{}'.format(val[0],val[1:])
    if len(val) == 4:
        t = '{}.{}'.format(val[0:2],val[2:])

    if j == 4:
        i,j = 1,0
        k = i+2

    for grp in list(fall.keys()):
        if grp  == 'PartType2': ### Disk for now
            type_all = fall['{}'.format(grp)]
            type_iso = fiso['{}'.format(grp)]
            posall = type_all['Coordinates']
            posiso = type_iso['Coordinates']

            sall_x,sall_y = posall[:,0],posall[:,1]
            siso_x,siso_y = posiso[:,0],posiso[:,1]

            axes[i,j].scatter(sall_x,sall_y,c='k',s=2,label='t={}'.format(t))
            axes[k,j].scatter(siso_x,siso_y,c='k',s=2,label='t={}'.format(t))
            axes[i,j].legend()
            axes[k,j].legend()

    j = j + 1

axes[0,0].set_ylabel('H+D+B')
axes[2,0].set_ylabel('D')
plt.subplots_adjust(hspace=0,wspace=0)
plt.show()
