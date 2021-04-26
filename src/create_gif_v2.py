### create a gif of v2 ###
### shamelessly stolen ###
### off stack exchange ### 
### differences only lie in
### location of plots & 
### location to place gif

import imageio
import glob
import argparse

parser = argparse.ArgumentParser(description =
    'Read Gadget Snapshots & plot')

parser.add_argument('combo',type = str,
    help = 'Halo,disk, bulge = all; disk= d') #v2 only allows disk or halo,disk
parser.add_argument('attemptN', type = int,   #v2 adds in an argument for attempt number
    help = 'Attempt Number')

args = parser.parse_args()

###v2 changes location 
path = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/previous_attempts/attempt{}/{}_results/plots_v2'.format(args.attemptN,args.combo)

all_files = glob.glob(path + '/*pos.png')
all_files2 = sorted(all_files)

images = []
for filename in all_files2:
    images.append(imageio.imread(filename))
### v2 changes location
imageio.mimsave('../previous_attempts/attempt{}/{}_results/plots_v2/test_combined.gif'.format(args.attemptN,args.combo),images,duration=0.3)


