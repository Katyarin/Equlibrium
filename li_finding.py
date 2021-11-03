import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 40805
time_list = []
W_list = []
betta_p_list = []
li_find_list = []
'''with open('W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))'''
W_list=[900]
time_list = [0.2]
start_time = t.time()
for ind in range(len(time_list)):
    print('______________________________________________________________________________________')
    print(ind)
    time = time_list[ind]
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/new_mcc/mcc_%d.json' % Shotn
    time, Rc = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    Globus_PET.COIL_upd(Shotn, time, I_coil)

    print(I_coil['Ipl'], Rc)
    betta_pol = 4 * W_list[ind] / (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rc)
    print('betta_pol: ', betta_pol)
    betta_p_list.append(betta_pol)

    betta_I = betta_pol / 1.2
    print('betta_I: ', betta_I)
    li = 1
    #alf11, alf22 = 1.5, 1
    li_f, res = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [60, 170, 10], False)
    li_prob = int(li_f * 100)
    li_acc, res_acc = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [li_prob - 10, li_prob + 10, 1], True)
    #print('real li: ', res['li'])
    #print('real bp: ', res['bp'])
    print(li_f, li_acc)
    betta_p_list.append(res_acc['bp'][res_acc['li'].index(li_acc)])
    li_find_list.append(li_acc)

with open('res_38800.txt', 'w') as res_file:
    for i in range(len(time_list)):
        res_file.write('%10.4f' % time_list[i])
        res_file.write('%10.4f' % W_list[i])
        res_file.write('%10.4f' % betta_p_list[i])
        res_file.write('%10.4f' % li_find_list[i])
        res_file.write('%10.4f' % (betta_p_list[i] + li_find_list[i] / 2))
        res_file.write('\n')

print(- start_time + t.time())







