import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
pi = 3.14159265359

Shotn = 38800
time = 0.18502
betta_po = 0.33
#betta_po = 0.9
alf11 = 1.18
alf22 = 0.965

dif_list = []
minimum = 1000
PATH = 'c:/work/equilibrium/Globus_PET/'
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

with open(PATH + 'shotn.txt', 'r') as file:
    for line in file:
        inf = line.split()
last_shot = int(inf[0])
time2 = float(inf[1])
Rc = float(inf[3])
#time = 0.18502
I_coil = {}
I_coil['Ipl'] = float(inf[2])
if last_shot != Shotn or time2 != time:
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/new_mcc/mcc_%d.json' % Shotn
    time, Rc = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    Globus_PET.COIL_upd(Shotn, time, I_coil)
    with open(PATH + 'shotn.txt', 'w') as file:
        file.write(' %s' %str(Shotn))
        file.write(' %s' %str(time))
        file.write(' %s' %str(I_coil['Ipl']))
        file.write(' %s' % str(Rc))

start_time = t.time()

#betta_find, res = Globus_PET.find_par('betta_po', Shotn, time, I_coil, betta_po, 1.27, [30, 50, 1])
betta_find = 0.249
'''print(res['bp'])
plt.figure()
plt.plot([i/100 for i in range(30, 30 + len(res['bp']))], [i/100 for i in range(30, 30 + len(res['bp']))])
plt.plot([i/100 for i in range(30, 30 + len(res['bp']))], res['bp'], 'r.')
plt.plot([i*1.2/100 for i in range(30, 30 + len(res['bp']))], res['bp'], 'g.')
plt.show()
print(betta_find)'''
li = 1
for g in range(5):
    li_find, res = Globus_PET.find_par('li', Shotn, time, I_coil, betta_find, li, [120, 121, 1])
print(t.time() - start_time)
print(res['li'])