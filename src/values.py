### Values that each file may use ###
### Keep & edit here              ###
### Execute file as needed        ###

M_disk = 1#5
M_halo = 1#12*M_disk
M_bulge = 1#M_disk/10
N_disk = 6000 #In Homogeneous modelI & HernI, it is 3000
N_halo = int(3*N_disk) #in homogeneous model I & hernI, it is 3*N_disk
N_bulge = int(N_disk/10) #all simulations
a_val = 1.0
h_val = 1.0
z_0 = 1.0
pb0,pd0,ph0 = 0.427,100.,0.711
mu_disk, sigma_disk = (1.0*10**(-3)), (1.0*10**(-5)) #HomoI, 1.0 & 0.1 for all
mu_halo, sigma_halo = (1.0*10**(-3)), (1.0*10**(-5))
mu_bulge, sigma_bulge = (1.0*10**(-3)), (1.0*10**(-5))


