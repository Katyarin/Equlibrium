import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t

pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 40269

time = 0.17
betta_I = 0.321
li = 0.79
I_coil_new = 288

try:
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
except FileNotFoundError:
    print('not found in new version')
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn

time, Rc, Rcm = Find_boundary.bound(f, time)
print('new time: ', time)

I_coil = Globus_PET.get_coils(Shotn, time)
print(I_coil)
if I_coil_new != 0:
    I_coil['Ipl'] = I_coil_new
print(I_coil)
Globus_PET.COIL_upd(Shotn, time, I_coil)
print(I_coil['Ipl'], Rc, Rcm)

Globus_PET.find_bound(Shotn, time, I_coil, betta_I, li, 0, 1, inside=True)
plt.show()