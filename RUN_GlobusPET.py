import Find_boundary
import numpy as np
import Globus_PET

Shotn = 38814
time = 0.187
betta_po = 0.46
alf11 = 1.04
alf22 = 1.1

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/mcc_%d.json' % Shotn

time = Find_boundary.bound(f, time)
print('new time: ', time)

PATH = 'c:/work/equilibrium/Globus_PET/'

I_coil = Globus_PET.get_coils(Shotn, time)
Globus_PET.DURS_upd(Shotn, time, I_coil['Ipl'], betta_po, alf11, alf22)
Globus_PET.COIL_upd(Shotn, time, I_coil)

