### test velocity dispersion function ###

import numpy as np
import functions as fcn

vArr = np.random.normal(10,1,100)
x = np.random.normal(1,0.1,100)
y = np.random.normal(2,0.2,100)
z = np.random.normal(3,0.3,100)

r = np.stack((x,y,z),axis=-1)
#print(r)
v_test,r_med_test = fcn.vel_disp(r,vArr,20)
#print(v_test,r_med_test)
