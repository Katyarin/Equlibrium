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

Shotn = 41103
mode = 'sep' #'sep' or 'psi'
dia_data = {}
try:
    p = 0
    with open('c:/work/equilibrium/dia_data/%d.txt' %Shotn, 'r') as file:
        for line in file:
            data = line.split(',')
            if p == 0:
                for i in data:
                    dia_data[i] = []
                p+=1
            else:
                for i, key in enumerate(list(dia_data.keys())):
                    if data[i] == '':
                        dia_data[key].append(data[i])
                    else:
                        dia_data[key].append(float(data[i]))
except:
    print('c:/work/equilibrium/dia_data/%d.txt' %Shotn)
    print('file not found')
    stop
print(dia_data)


I_coil_new = dia_data['Ip']
time_list = []
W_list = []
betta_p_list = []
betta_I_list = []
li_find_list = []
li_code_list = []
W_all_list = []
We_list = []
Wi_list = []
V_list = []
S_list = []
p_axis_list = []
ftorpl_list = []

pdf_file = PdfPages(str(Shotn) + '_plots.pdf')

'''with open('files/W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))
'''

#W_list=[8600]
#betta_I_list_raw = [0.28, 0.28, 0.28, 0.29, 0.31, 0.33, 0.36, 0.39, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.46, 0.47, 0.43, 0.39]
time_list = [i/1000 for i in dia_data['time']]
betta_I_list = [round(i, 2) for i in dia_data['betadia\n']]
Bt = dia_data['Bt']
print(I_coil_new, time_list, betta_I_list, Bt)
print(len(I_coil_new), len(time_list),len(betta_I_list),len(Bt))

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
    #f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    time, Rc, Rcm = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    print(I_coil)
    if I_coil_new:
        if I_coil_new[ind] != 0:
            I_coil['Ipl'] = I_coil_new[ind] * 1e3
    print(I_coil)

    Globus_PET.COIL_upd(Shotn, time, I_coil, Bt[ind])

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
    li_f, li_f2, res, psi_delta, bounds_delta = Globus_PET.find_par2('li', mode, Shotn, time, I_coil, betta_I, li, [50, 180, 10], pdf_file, True, 'all', psi_dia_sakh=dia_data['psidia'][ind])
    if mode == 'sep':
        li_prob = int(li_f2 * 100)
        li_acc, li_acc2, res_acc, psi_delta_acc, bounds_delta_acc  = Globus_PET.find_par2('li', mode,  Shotn, time, I_coil, betta_I, li, [li_prob - 7, li_prob + 7, 1], pdf_file, True, 'all')
    elif mode == 'psi':
        li_acc2 = li_f2
    else:
        print('unknown mode!')
        stop
    b_I, li_3, bp, W_all, We, V, S, P_axis, li_code, Ftor_pl, Wi = Globus_PET.find_bound(Shotn, time, I_coil, betta_I, li_acc2, 0, 1, pdf = pdf_file,
                                                                        inside=True, Wi=True, share=True)
    #plt.show()
    #print('real li: ', res['li'])
    #print('real bp: ', res['bp'])
    #print(li_f, li_acc)
    betta_p_list.append(float(bp))
    li_find_list.append(float(li_3))
    li_code_list.append(float(li_code))
    W_all_list.append(W_all)
    We_list.append(We)
    V_list.append(V)
    S_list.append(S)
    p_axis_list.append(P_axis)
    ftorpl_list.append(Ftor_pl)
    Wi_list.append(Wi)
    #show = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [int(li_acc2 * 100), int(li_acc2 * 100) + 1, 1], True, 'all')

pdf_file.close()

with open('220722_6res_' + str(Shotn) +'.txt', 'w') as res_file:
    for i in range(delta_index):
        res_file.write('%8.4f'% time_list[i])
        res_file.write('%8.4f'% betta_I_list[i])
        res_file.write('%8.4f'% betta_p_list[i])
        res_file.write('%8.4f'% li_find_list[i])
        res_file.write('%14.4f' % W_all_list[i])
        res_file.write('%14.4f' % W_list[i])
        res_file.write('%14.4f' % We_list[i])
        res_file.write('%8.4f'% V_list[i])
        res_file.write('%8.4f'% S_list[i])
        res_file.write('%8.4f'% p_axis_list[i])
        res_file.write('%8.4f'% ftorpl_list[i])
        res_file.write('%8.4f' % li_code_list[i])
        res_file.write('%14.4f' % Wi_list[i])
        res_file.write('\n')

#plt.show()

print('time left: ', - start_time + t.time())







