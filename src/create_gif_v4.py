### create a gif of v3 ###
### shamelessly stolen ###
### had to downsize gifs 
### because files wouldn't fit

import imageio
import glob

aN = 16
combo = 'all'
foldername = 'plots_v3'
fname = 'hern_approx_r_all_v3.gif'

path = r'/N/u/bkkimsey/Carbonate/gadget/Gadget-2.0.7/Gadget2/previous_attempts/attempt{}/{}_results/{}'.format(aN,combo,foldername)

all_files = glob.glob(path + '/*pos.png')
all_files2 = sorted(all_files)

images = []

for i,filename in enumerate(all_files2):
    if i<=150:
        if (i % 2) == 0:
            images.append(imageio.imread(filename))
     
imageio.mimsave('../previous_attempts/attempt{}/{}_results/{}/{}'.format(aN,combo,foldername,fname),images,duration=0.3)

