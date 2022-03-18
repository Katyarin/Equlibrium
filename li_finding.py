import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
from matplotlib.backends.backend_pdf import PdfPages



pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 40269
I_coil_new = []
time_list = []
W_list = []
betta_p_list = []
betta_I_list = []
li_find_list = []
li_code_list = []

pdf_file = PdfPages(str(Shotn) + '_plots.pdf')

'''with open('files/W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))
'''

#W_list=[8600]
betta_I_list = [0.321]
time_list = [0.17]


index_start = 0
#index_start = 23
index_end = len(time_list)
#index_end = 24

delta_index = index_end - index_start
start_time = t.time()
for ind in range(index_start, index_end):
    print('______________________________________________________________________________________')
    print(ind)
    time = time_list[ind]
    try:
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
    except FileNotFoundError:
        print('not found in new version')
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    '''if Shotn < 41033:
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    else:
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn'''
    time, Rc, Rcm = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    print(I_coil)
    if I_coil_new:
        if I_coil_new[ind] != 0:
            I_coil[ind]['Ipl'] = I_coil_new
    print(I_coil)

    Globus_PET.COIL_upd(Shotn, time, I_coil)

    print(I_coil['Ipl'], Rc, Rcm)
    '''betta_pol = 4 * W_list[ind] / (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rc)
    print('betta_pol: ', betta_pol)
    betta_p_list.append(betta_pol)'''

    if len(W_list) > ind:
        betta_I = 4 * W_list[ind] / (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rcm)
        print('betta_I: ', betta_I)
        betta_I_list.append(betta_I)
    elif len(betta_I_list) > ind:
        betta_I = betta_I_list[ind]
        W = betta_I * (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rcm) / 4
        W_list.append(W)
        print('W = ', W)
    else:
        print('SOMETHING WRONG')
        stop
    li = 1
    #alf11, alf22 = 1.5, 1
    li_f, li_f2, res = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [50, 180, 10], pdf_file, True, 'all')
    li_prob = int(li_f2 * 100)
    li_acc, li_acc2, res_acc = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [li_prob - 7, li_prob + 7, 1], pdf_file, True, 'all')
    #plt.show()
    #print('real li: ', res['li'])
    #print('real bp: ', res['bp'])
    print(li_f, li_acc)
    betta_p_list.append(res_acc['bp'][res_acc['li'].index(li_acc)])
    li_find_list.append(li_acc)
    li_code_list.append(li_acc2)
    #show = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [int(li_acc2 * 100), int(li_acc2 * 100) + 1, 1], True, 'all')

pdf_file.close()

with open('res_' + str(Shotn) +'.txt', 'w') as res_file:
    for i in range(delta_index):
        res_file.write('%10.4f' % time_list[i + index_start])
        res_file.write('%14.4f' % W_list[i + index_start])
        res_file.write('%10.4f' % betta_I_list[i])
        res_file.write('%10.4f' % betta_p_list[i])
        res_file.write('%10.4f' % li_find_list[i])
        res_file.write('%10.4f' % (betta_p_list[i] + li_find_list[i] / 2))
        res_file.write('\n')

#plt.show()

print(- start_time + t.time())







