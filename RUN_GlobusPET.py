import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
from matplotlib.backends.backend_pdf import PdfPages
import json
import datetime

b = str(datetime.date.today())
tomorow = b[2:4] + b[5:7] + b[8:]
#tomorow = '230329'
date_of_culc = '240604'
Path_res = 'results/'

shot_list = [44330]
pf2=1
for Shotn in shot_list:
    print(Shotn)
    shot_without_plasma = 0
    mode = 'psi' #'sep' or 'psi'
    my = 1
    coincidence = True

    pi = 3.14159265359
    mu0 = 4*pi*1e-7
    PATH = 'c:/work/equilibrium/Globus_PET/'
    path_res = 'eq_results/%s/' %Shotn

    with open('eq_results/statistic.json', 'r') as statFile:
        statistic = json.load(statFile)
    if str(Shotn) in statistic['data_of_culc'].keys() and date_of_culc == False:
        date_of_culc = statistic['data_of_culc'][str(Shotn)]
        print('we here?')
    else:
        statistic['data_of_culc'][Shotn] = date_of_culc
        print('or here?')
    print(date_of_culc)
    if shot_without_plasma:
        statistic['shot_without_plasma'][Shotn] = shot_without_plasma

    stat_list ={}
    for key in statistic:
        if key != 'data_of_culc' and key != 'shot_without_plasma':
            statistic[key][Shotn] = {'av': 0, 'range': 0}
            stat_list[key] = []

    dia_data = Globus_PET.openDiaFile(Shotn, my)

    print(dia_data)
    I_coil_new = dia_data['Ip']
    print(len(I_coil_new))
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

    alpha = []
    if 'alpha1' in culc_data.keys():
        alpha = [[culc_data['alpha1'][i], culc_data['alpha2'][i]] for i in range(len(culc_data['alpha1']))]
        print('yes alpha!!')
    else:
        print('no alpha =(')

    with open(path_res + tomorow + '_res_' + str(Shotn) + '.txt', 'w') as res_file:
        res_file.write('%8s' % 'time')
        res_file.write('%8s' % 'beta_I')
        res_file.write('%8s' % 'beta_p')
        res_file.write('%8s' % 'li')
        res_file.write('%14s' % 'W')
        res_file.write('%14s' % 'W_approxRav')
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
        res_file.write('%14s' % 'beta_t')
        res_file.write('\n')

    print(len(time_list))

    for ind, time in enumerate(time_list):
        try:
            f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
        except FileNotFoundError:
            print('not found in new version')
        # f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
        time, Rc, Rcm, Rav, k, conf = Find_boundary.bound(f, time)
        conf='lim'
        print('new time: ', time)

        I_coil = Globus_PET.get_coils(Shotn, time)

        print(I_coil)
        if I_coil_new:
            if I_coil_new[ind] != 0:
                I_coil['Ipl'] = I_coil_new[ind] * 1e3

        '''plasma shift!!!!!!!!!!!!!!!!!!!!!!
        I_coil['Ipf2'] = 1000 *1e-3'''
        print(I_coil)

        betta_I = betta_I_list[ind]
        li_acc2 = li_list[ind]
        if alpha:
            alpha_loc = alpha[ind]
        else:
            alpha_loc = li_acc2

        #Globus_PET.COIL_upd(Shotn, time, I_coil, Bt[ind], Rav, Ipf2oposite=True)
        if pf2:
            coil_data = {
                'PF1_0': {'Icoil': I_coil['Ipf1'] * 1e-3},
                'PF1_1': {'Icoil': I_coil['Ipf1'] * 1e-3},
                'PF2_0': {'Icoil': dia_data['pf2Up'][ind] * 1e-3},
                'PF2_1': {'Icoil': I_coil['Ipf2'] * 1e-3, 'link': 3},
                'PF3_0': {'Icoil': I_coil['Ipf3'] * 1e-3},
                'PF3_1': {'Icoil': I_coil['Ipf3'] * 1e-3},
                'HFC1_0': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC1_1': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC2_0': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC2_1': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'VFC_0': {'Icoil': I_coil['Ivfc'] * 1e-3},
                'VFC_1': {'Icoil': I_coil['Ivfc'] * 1e-3},
                'CC1_0': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 7},
                'CC1_1': {'Icoil': I_coil['Icc'] * 1e-3},
                'CC2_0': {'Icoil': I_coil['Icc'] * 1e-3},
                'CC2_1': {'Icoil': I_coil['Icc'] * 1e-3},
                'CC3_0': {'Icoil': I_coil['Icc'] * 1e-3},
                'CC3_1': {'Icoil': I_coil['Icc'] * 1e-3},
                'CS_0': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8},
                'CS_1': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8},
                'CS_2': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8},
                'CS_3': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8},
                'CS_4': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8},
                'CS_5': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 8}
            }
            # print(fast_ind)
            Globus_PET.newCOIL_upd(Shotn, time, I_coil['Ipl'], Rav * 100, Bt[ind], coil_data)
        else:
            coil_data = {
                'PF1_0': {'Icoil': I_coil['Ipf1'] * 1e-3},
                'PF1_1': {'Icoil': I_coil['Ipf1'] * 1e-3},
                'PF2_0': {'Icoil': I_coil['Ipf2'] * 1e-3},
                'PF2_1': {'Icoil': I_coil['Ipf2'] * 1e-3, 'link': 2},
                'PF3_0': {'Icoil': I_coil['Ipf3'] * 1e-3},
                'PF3_1': {'Icoil': I_coil['Ipf3'] * 1e-3},
                'HFC1_0': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC1_1': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC2_0': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'HFC2_1': {'Icoil': I_coil['Ihfc'] * 1e-3},
                'VFC_0': {'Icoil': I_coil['Ivfc'] * 1e-3},
                'VFC_1': {'Icoil': I_coil['Ivfc'] * 1e-3},
                'CC1_0': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CC1_1': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CC2_0': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CC2_1': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CC3_0': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CC3_1': {'Icoil': I_coil['Icc'] * 1e-3, 'link': 6},
                'CS_0': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7},
                'CS_1': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7},
                'CS_2': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7},
                'CS_3': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7},
                'CS_4': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7},
                'CS_5': {'Icoil': I_coil['Ics'] * 1e-3, 'link': 7}
            }
            # print(fast_ind)
            Globus_PET.newCOIL_upd(Shotn, time, I_coil['Ipl'], Rav * 100, Bt[ind], coil_data)
        Globus_PET.DATA_upd(conf)
        bounds_delta_res = 1000
        trying = 0
        coid = True
        if 'boundDr' not in culc_data:
            culc_data['boundDr'] = [1]*len(time_list)
        while abs(float(bounds_delta_res) - float(culc_data['boundDr'][ind])) > 0.00005 and coid:
            b_I, li_3, bp, W_all, We, V, S, P_axis, li_code, Ftor_pl, Wi, bounds_delta_res, ne_av, strike_point, q95 = Globus_PET.find_bound(Shotn,
                                                                                                                      time,
                                                                                                                      I_coil,
                                                                                                                      betta_I,
                                                                                                                      alpha_loc,
                                                                                                                      0, 1,
                                                                                                                      pdf=pdf_file,
                                                                                                                      inside=True,
                                                                                                                      Wi=True,
                                                                                                                      share=True)
            coid = coincidence
            trying +=1
            print(li_acc2)
            if trying < 6:
                li_acc2 += 0.00001
            elif trying == 6:
                li_acc2 = li_acc2 - 0.00006
            else:
                li_acc2 = li_acc2 - 0.00001
            print(li_acc2)
            print('time: ', time, ', trying no ', trying, ' ', abs(float(bounds_delta_res) - float(culc_data['boundDr'][ind])))
        rb = []
        zb = []
        with open(Path_res + str(Shotn) + '_' + str(round(time, 3)) + '_bound.txt', 'r') as file2:
            for line in file2:
                data = line.split()
                rb.append(float(data[0]))
                zb.append((float(data[1])))
        z_up = max(zb)
        z_low = min(zb)
        print(z_low)
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
        print(dots['x-dot'][1])
        if z_low < dots['x-dot'][1]:
            if z_low*dots['x-dot'][1]>0:
                z_low = dots['x-dot'][1]
                r_low = dots['x-dot'][0]

        R = (r_out + r_in)/2
        Rav = Rav/100
        a = (r_out - r_in)/2
        b = (z_up - z_low) /2
        print(r_out, r_in, z_up, z_low)
        print(a, b)
        k = b/a
        print(k)
        tr_up = (R-r_low) /a
        tr_down = (R-r_up) /a
        stat_list['beta_I'].append(float(betta_I))
        stat_list['li'].append(float(li_3))
        stat_list['Bt'].append(Bt[ind])
        stat_list['Ip'].append(I_coil['Ipl'])
        stat_list['W'].append(W_all)
        stat_list['We'].append(We)
        if Wi != 0:
            stat_list['Wi'].append(Wi)
        stat_list['<ne>'].append(ne_av)
        stat_list['k'].append(k)
        stat_list['A'].append(R/a)
        with open(path_res + tomorow + '_res_' + str(Shotn) + '.txt', 'a') as res_file:
            res_file.write('%8.4f' % time)
            res_file.write('%8.4f' % float(betta_I))
            res_file.write('%8.4f' % float(bp))
            res_file.write('%8.4f' % float(li_3))
            res_file.write('%14.4f' % W_all)
            res_file.write('%14.4f' % (3/2 *betta_I_sakharov[ind] * (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rav) / 4))
            res_file.write('%14.4f' % We)
            res_file.write('%8.4f' % V)
            res_file.write('%8.4f' % S)
            res_file.write('%8.4f' % P_axis)
            res_file.write('%8.5f' % Ftor_pl)
            res_file.write('%8.4f' % float(li_acc2))
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
            res_file.write('%14.4f' % (1.6 * pi * W_all *1e-4 / (3*Bt[ind]*Bt[ind]*V)))
            res_file.write('\n')
        #cont = int(input('continue?'))
    pdf_file.close()
    for key in stat_list.keys():
        print(key)
        if stat_list[key]:
            statistic[key][Shotn]['av'] = sum(stat_list[key])/len(stat_list[key])
            statistic[key][Shotn]['range'] = (max(stat_list[key]) - min(stat_list[key]))/2
        else:
            print('no data')

    with open('eq_results/statistic.json', 'w') as statFileres:
        json.dump(statistic, statFileres)