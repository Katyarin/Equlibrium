import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
from matplotlib.backends.backend_pdf import PdfPages

tomorow = '221202'
date_of_culc = '221202_2'
Path_res = 'results/'

Shotn = 42567
mode = 'psi' #'sep' or 'psi'

pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'
path_res = 'eq_results/%s/' %Shotn

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

print(dia_data)
I_coil_new = dia_data['Ip']
time_list = [i / 1000 for i in dia_data['time']]
betta_I_sakharov = [round(i, 2) for i in dia_data['betadia\n']]
betta_I_list = []
Bt = dia_data['Bt']

culc_data = {}
data_name = []
l1 = 0
with open(path_res + 'test_' + date_of_culc + 'res_' + str(Shotn) + '.txt', 'r') as res_file:
    for line in res_file:
        if l1 == 0:
            for data_key in line.split():
                culc_data[data_key] = []
                data_name.append(data_key)
        else:
            for i, j in enumerate(line.split()):
                culc_data[data_name[i]].append(float(j))
        l1 +=1
print(culc_data)
time_list = culc_data['time']
betta_I_list = culc_data['beta_I']
li_list = culc_data['li_code']
pdf_file = PdfPages(path_res + tomorow + '_' + str(Shotn) + '_second_plots.pdf')

with open(path_res + 'test_' + tomorow + 'res_2_' + str(Shotn) + '.txt', 'a') as res_file:
    res_file.write('%8s' % 'time')
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
    res_file.write('%8s' % 'boundDr')
    res_file.write('%8s' % 'Ip')
    res_file.write('%8s' % 'Bt')
    res_file.write('%14s' % '<ne>')
    res_file.write('%8s' % 'R')
    res_file.write('%8s' % 'a')
    res_file.write('%8s' % 'r_ax')
    res_file.write('%8s' % 'k')
    res_file.write('%8s' % 'tr_up')
    res_file.write('%8s' % 'tr_down')
    res_file.write('%14s' % 'r_sp_in')
    res_file.write('%14s' % 'z_sp_in')
    res_file.write('%14s' % 'r_sp_out')
    res_file.write('%14s' % 'z_sp_out')
    res_file.write('%14s' % 'q95')
    res_file.write('\n')

for ind, time in enumerate(time_list):
    try:
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
    except FileNotFoundError:
        print('not found in new version')
    # f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    time, Rc, Rcm = Find_boundary.bound(f, time)
    print('new time: ', time)

    I_coil = Globus_PET.get_coils(Shotn, time)
    print(I_coil)
    if I_coil_new:
        if I_coil_new[ind] != 0:
            I_coil['Ipl'] = I_coil_new[ind] * 1e3

    betta_I = betta_I_list[ind]
    li_acc2 = li_list[ind]

    Globus_PET.COIL_upd(Shotn, time, I_coil, Bt[ind])
    b_I, li_3, bp, W_all, We, V, S, P_axis, li_code, Ftor_pl, Wi, bounds_delta_res, ne_av, strike_point, q95 = Globus_PET.find_bound(Shotn,
                                                                                                                  time,
                                                                                                                  I_coil,
                                                                                                                  betta_I,
                                                                                                                  li_acc2,
                                                                                                                  0, 1,
                                                                                                                  pdf=pdf_file,
                                                                                                                  inside=True,
                                                                                                                  Wi=True,
                                                                                                                  share=True)
    rb = []
    zb = []
    with open(Path_res + str(Shotn) + '_' + str(round(time, 3)) + '_bound.txt', 'r') as file2:
        for line in file2:
            data = line.split()
            rb.append(float(data[0]))
            zb.append((float(data[1])))
    z_up = max(zb)
    z_low = min(zb)
    r_out = max(rb)
    r_in = min(rb)
    r_up = rb[zb.index(z_up)]
    r_low = rb[zb.index(z_low)]
    dots = {}
    with open(Path_res + str(Shotn) + '_' + str(round(time, 3)) + '_dots.txt', 'r') as file3:
        for line in file3:
            data = line.split()
            if len(data) < 4:
                dots[data[0]] = [float(data[1]), float(data[2])]
            else:
                dots[data[0] + '_' + data[1]] = [float(data[2]), float(data[3])]
    r_ax = dots['magnetic_axis'][0]
    if z_low < dots['x-dot'][1]:
        z_low = dots['x-dot'][1]
        r_low = dots['x-dot'][0]

    R = (r_out+ r_in)/2
    a = (r_out- r_in)/2
    b = (z_up - z_low) /2
    k = b/a
    tr_up = (R-r_low) /a
    tr_down = (R-r_up) /a
    with open(path_res + 'test_' + tomorow + 'res_2_' + str(Shotn) + '.txt', 'a') as res_file:
        res_file.write('%8.4f' % time)
        res_file.write('%8.4f' % float(betta_I))
        res_file.write('%8.4f' % float(bp))
        res_file.write('%8.4f' % float(li_3))
        res_file.write('%14.4f' % W_all)
        res_file.write('%14.4f' % (3/2 *betta_I * (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*r_ax) / 4))
        res_file.write('%14.4f' % We)
        res_file.write('%8.4f' % V)
        res_file.write('%8.4f' % S)
        res_file.write('%8.4f' % P_axis)
        res_file.write('%8.5f' % Ftor_pl)
        res_file.write('%8.4f' % float(li_code))
        res_file.write('%14.4f' % Wi)
        res_file.write('%8.4f' % float(bounds_delta_res))
        res_file.write('%8.2f' % I_coil['Ipl'])
        res_file.write('%8.4f' % Bt[ind])
        res_file.write('%14.4e' % ne_av)
        res_file.write('%8.4f' % R)
        res_file.write('%8.4f' % a)
        res_file.write('%8.4f' % r_ax)
        res_file.write('%8.4f' % k)
        res_file.write('%8.4f' % tr_up)
        res_file.write('%8.4f' % tr_down)
        res_file.write('%14.4f' % strike_point['inner'][0])
        res_file.write('%14.4f' % strike_point['inner'][1])
        res_file.write('%14.4f' % strike_point['outer'][0])
        res_file.write('%14.4f' % strike_point['outer'][1])
        res_file.write('%14.4f' % q95)
        res_file.write('\n')

pdf_file.close()