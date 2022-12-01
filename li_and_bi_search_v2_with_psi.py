import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
import json
from matplotlib.backends.backend_pdf import PdfPages
import os

#this test for search equlibrium with beta and li

tomorow = '221201'

Shotn = 42567
mode_start = 'psi' #'sep' or 'psi'
current_by_FreeGS = False

pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'
path_res = 'eq_results/%s/' %Shotn

try:
    os.mkdir(path_res)
except OSError:
    print('Не удалось создать папку')

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
                    if data[i]:
                        dia_data[key].append(float(data[i]))
                    else:
                        dia_data[key].append(0)
except:
    print('c:/work/equilibrium/dia_data/%d.txt' %Shotn)
    print('file not found')
    stop
print(dia_data)

try:
    l = 0
    with open('c:/work/equilibrium/Ti_data/' + str(Shotn) + '.txt', 'r') as Ti_file:
        for line in Ti_file:
            data = line.split()
            if l == 0:
                time_li = [float(data[i]) for i in range(1, len(data), 2)]
            l+=1
    print(time_li)
    time_min = min(dia_data['time'])
    time_max = max(dia_data['time'])
    for t1 in time_li:
        if time_min < t1 < time_max:
            if t1 not in dia_data['time']:
                for i,t2 in enumerate(dia_data['time']):
                    if t1 < t2:
                        index = i
                        break
                for key in dia_data.keys():
                    if key=='time':
                        dia_data[key].insert(index, t1)
                    else:
                        t2 = dia_data['time'][index+1]
                        t0 = dia_data['time'][index-1]
                        v2 = dia_data[key][index]
                        v0 = dia_data[key][index-1]
                        value = v0 + (v2-v0)*(t1-t0) / (t2 - t0)
                        dia_data[key].insert(index, value)
except FileNotFoundError:
    print("NO Ti data!")
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

pdf_file =PdfPages(path_res + tomorow + '_' +str(Shotn) + '_plots.pdf')

'''with open('files/W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))
'''

#W_list=[8600]
#betta_I_list_raw = [0.28, 0.28, 0.28, 0.29, 0.31, 0.33, 0.36, 0.39, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.46, 0.47, 0.43, 0.39]
time_list = [i/1000 for i in dia_data['time']]
betta_I_sakharov = [round(i, 2) for i in dia_data['betadia\n']]
betta_I_list = []
Bt = dia_data['Bt']
print(I_coil_new, time_list, betta_I_list, Bt)
print(len(I_coil_new), len(time_list),len(betta_I_list),len(Bt))

index_start = 0
#index_start = 23
index_end = len(time_list)
#index_end = 24

delta_index = index_end - index_start
start_time = t.time()
time_already_done = []
try:
    l1 = 0
    with open(path_res + 'test_' + tomorow + 'res_' + str(Shotn) + '.txt', 'r') as res_file:
        for line in res_file:
            if l1 != 0:
                time_already_done.append(float(line.split()[0]))

            l1 +=1
except FileNotFoundError:
    with open(path_res + 'test_' + tomorow + 'res_' + str(Shotn) +'.txt', 'a') as res_file:
        res_file.write('%8s' %'time')
        res_file.write('%8s' % 'beta_I')
        res_file.write('%8s' % 'beta_p')
        res_file.write('%8s' % 'li')
        res_file.write('%14s' % 'W')
        res_file.write('%14s' % 'W_approx')
        res_file.write('%14s' % 'We')
        res_file.write('%8s' % 'V')
        res_file.write('%8s' % 'S')
        res_file.write('%8s' % 'P_axis')
        res_file.write('%8s' % 'Psi')
        res_file.write('%8s' % 'li_code')
        res_file.write('%14s' % 'Wi')
        res_file.write('%8s' % 'boundD')
        res_file.write('%8s' % 'boundDr')
        res_file.write('%8s' % 'Ip')
        res_file.write('%8s' % 'Bt')
        res_file.write('%14s' % '<ne>')
        res_file.write('\n')

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

    if round(time,4) in time_already_done:
        print('this time alredy done')
        continue

    I_coil = Globus_PET.get_coils(Shotn, time)
    print(I_coil)
    if I_coil_new:
        if I_coil_new[ind] != 0:
            I_coil['Ipl'] = I_coil_new[ind] * 1e3
    #print(I_coil)

    if current_by_FreeGS:
        with open('currents_by_FreeGS/wall_and_loop_'+ str(Shotn) + '_00' + str(int(time*1000)) + '.json', 'r') as file_GS:
            data_zhenya = json.load(file_GS)
        FreeGS_curr = data_zhenya['coil_signals']

        for curr in I_coil.keys():
            I_coil[curr] = FreeGS_curr[curr + '_sm'] / 1000
        #print(I_coil)

    Globus_PET.COIL_upd(Shotn, time, I_coil, Bt[ind])
    #print(I_coil['Ipl'], Rc, Rcm)

    li_for_one_time = []
    li_work_one_time = []
    bI_for_one_time = []
    psi_delta_for_one_time = []
    bounds_delta_for_one_time = []
    bound_r = (int(betta_I_sakharov[ind]*100) + 10)
    bound_l = (int(betta_I_sakharov[ind]*100) - 10)
    mode = mode_start
    beta_0 = betta_I_sakharov[ind] - 0.05
    psi_delta = 1
    bounds_delta = 1
    attempt = 0
    while (psi_delta**2 + bounds_delta**2)**0.5 > 0.0065:
        if attempt % 2 and attempt < 10:
            betta_I = beta_0 - (attempt%2 + attempt//2) * 0.01
        else:
            betta_I = beta_0 + (attempt % 2 + attempt // 2) * 0.01
        print('attempt #%i, beta=%f' %(attempt, betta_I))
        bI_for_one_time.append(betta_I)
        li = 1
        #alf11, alf22 = 1.5, 1
        li_f, li_f2, res, psi_delta, bounds_delta = Globus_PET.find_par2('li', mode, Shotn, time, I_coil, betta_I, li, [50, 180, 10], pdf_file, False, 'all', psi_dia_sakh=dia_data['psidia'][ind])

        li_for_one_time.append(li_f)
        li_work_one_time.append(li_f2)
        psi_delta_for_one_time.append(psi_delta)
        bounds_delta_for_one_time.append(bounds_delta)
        if attempt<10:
            attempt+=1
        elif attempt<30:
            attempt+=2
        else:
            break
    '''pdf_file.close()
    plt.figure()
    plt.plot(bI_for_one_time, psi_delta_for_one_time)
    plt.grid()
    plt.ylabel('delta_psi')
    plt.figure()
    plt.plot(bI_for_one_time, li_for_one_time)
    plt.grid()
    plt.ylabel('li')
    plt.figure()
    plt.plot(bI_for_one_time, bounds_delta_for_one_time)
    plt.grid()
    plt.ylabel('delta_bounds')
    plt.show()'''
    ind_min = bounds_delta_for_one_time.index(min(bounds_delta_for_one_time))
    print('before check')
    print(bounds_delta_for_one_time[ind_min], psi_delta_for_one_time[ind_min])
    mera_list1 = []
    mera_list2 = []
    if psi_delta_for_one_time[ind_min] > 0.01:
        for i in range(len(bounds_delta_for_one_time)):
            mera_list1.append((psi_delta_for_one_time[i] + bounds_delta_for_one_time[i])/2)
            mera_list2.append((psi_delta_for_one_time[i] ** 2 + bounds_delta_for_one_time[i] **2) **0.5)
        print(min(mera_list1), mera_list1.index(min(mera_list1)),bounds_delta_for_one_time[mera_list1.index(min(mera_list1))], psi_delta_for_one_time[mera_list1.index(min(mera_list1))])
        print(min(mera_list2), mera_list2.index(min(mera_list2)),bounds_delta_for_one_time[mera_list2.index(min(mera_list2))],
              psi_delta_for_one_time[mera_list2.index(min(mera_list2))])
        ind_min = mera_list1.index(min(mera_list1))
    print('after check')
    print(bounds_delta_for_one_time[ind_min], psi_delta_for_one_time[ind_min])
    betta_I = bI_for_one_time[ind_min]
    li_acc2 = li_work_one_time[ind_min]
    b_I, li_3, bp, W_all, We, V, S, P_axis, li_code, Ftor_pl, Wi, bounds_delta_res, ne_av, strike_point = Globus_PET.find_bound(Shotn, time, I_coil, betta_I, li_acc2, 0, 1, pdf = pdf_file,
                                                                        inside=True, Wi=True, share=True)
    #plt.show()
    #print('real li: ', res['li'])
    #print('real bp: ', res['bp'])
    #print(li_f, li_acc)
    betta_p_list.append(float(bp))
    betta_I_list.append(float(betta_I))
    li_find_list.append(float(li_3))
    li_code_list.append(float(li_code))
    W_all_list.append(W_all)
    We_list.append(We)
    V_list.append(V)
    S_list.append(S)
    p_axis_list.append(P_axis)
    ftorpl_list.append(Ftor_pl)
    Wi_list.append(Wi)
    W_list.append(betta_I * (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rcm) / 4)
    with open(path_res + 'test_' + tomorow + 'res_' + str(Shotn) + '.txt', 'a') as res_file:
        res_file.write('%8.4f' % time)
        res_file.write('%8.4f' % float(betta_I))
        res_file.write('%8.4f' % float(bp))
        res_file.write('%8.4f' % float(li_3))
        res_file.write('%14.4f' % W_all)
        res_file.write('%14.4f' % (betta_I * (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rcm) / 4))
        res_file.write('%14.4f' % We)
        res_file.write('%8.4f' % V)
        res_file.write('%8.4f' % S)
        res_file.write('%8.4f' % P_axis)
        res_file.write('%8.5f' % Ftor_pl)
        res_file.write('%8.4f' % float(li_code))
        res_file.write('%14.4f' % Wi)
        res_file.write('%8.4f' % float(bounds_delta_for_one_time[ind_min]))
        res_file.write('%8.4f' % float(bounds_delta_res))
        res_file.write('%8.2f' % I_coil['Ipl'])
        res_file.write('%8.4f' % Bt[ind])
        res_file.write('%14.4e' % ne_av)
        res_file.write('\n')
    #show = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [int(li_acc2 * 100), int(li_acc2 * 100) + 1, 1], True, 'all')


pdf_file.close()

with open(path_res + tomorow + 'res_' + str(Shotn) +'.txt', 'w') as res_file:
    res_file.write('%8s' %'time')
    res_file.write('%8s' % 'beta_I')
    res_file.write('%8s' % 'beta_p')
    res_file.write('%8s' % 'li')
    res_file.write('%14s' % 'W')
    res_file.write('%14s' % 'W_approx')
    res_file.write('%14s' % 'We')
    res_file.write('%8s' % 'V')
    res_file.write('%8s' % 'S')
    res_file.write('%8s' % 'P_axis')
    res_file.write('%8s' % 'Psi')
    res_file.write('%8s' % 'li_code')
    res_file.write('%14s' % 'Wi')
    res_file.write('\n')
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
        res_file.write('%8.5f'% ftorpl_list[i])
        res_file.write('%8.4f' % li_code_list[i])
        res_file.write('%14.4f' % Wi_list[i])
        res_file.write('\n')

#plt.show()

print('time left: ', - start_time + t.time())




