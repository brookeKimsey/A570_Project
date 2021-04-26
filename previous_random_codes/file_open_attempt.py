### first attempt to load in .hdf5 files

import h5py

f = h5py.File('mytestfile.hdf5','r')

print(list(f.keys()))

dset1 = f['Type1/Positions']
print(dset1[0])
