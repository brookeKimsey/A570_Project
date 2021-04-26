###    create a gif of original plots    ###
### shamelessly stolen off stackexchange ### 

import imageio
import glob
import argparse

parser = argparse.ArgumentParser(description =
    'Read Gadget Snapshots & plot')

parser.add_argument('combo',type = str,
    help = 'Halo,disk, bulge = all; halo,disk=hd; disk,bulge = bd; disk= d')

args = parser.parse_args()

path = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/attempt/{}_results/plots'.format(args.combo)

all_files = glob.glob(path + '/*pos.png')
all_files2 = sorted(all_files)

images = []
for filename in all_files2:
    images.append(imageio.imread(filename))
imageio.mimsave('../attempt/{}_results/plots/test_combined.gif'.format(args.combo),images,duration=0.3)


