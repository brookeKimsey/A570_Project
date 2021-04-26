### create a gif of v3 ###
### shamelessly stolen ###
### changed from original
### version to get images
### and put gif from/in plots_v3  

import imageio
import glob
import argparse

parser = argparse.ArgumentParser(description =
    'Read Gadget Snapshots & plot')

parser.add_argument('combo',type = str,
    help = 'Halo,disk, bulge = all; halo,disk=hd; disk,bulge = bd; disk= d')
parser.add_argument('attemptN', type = int,   #v2 adds in an argument for attempt number
    help = 'Attempt Number')

args = parser.parse_args()

path = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/previous_attempts/attempt{}/{}_results/plots_v3'.format(args.attemptN,args.combo)

all_files = glob.glob(path + '/*pos.png')
all_files2 = sorted(all_files)

images = []
for filename in all_files2:
    images.append(imageio.imread(filename))
imageio.mimsave('../previous_attempts/attempt{}/{}_results/plots_v3/test_combined.gif'.format(args.attemptN,args.combo),images,duration=0.3)


