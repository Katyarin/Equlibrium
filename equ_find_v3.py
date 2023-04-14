import Find_boundary
import Globus_PET
import matplotlib.pyplot as plt
import time as t
import json
from matplotlib.backends.backend_pdf import PdfPages
import os
import get_dia_data
import datetime
pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

#this test for fast search equlibrium with beta and li by delta_bound and delta_psi

my = 1
b = str(datetime.date.today())
tomorow = b[2:4] + b[5:7] + b[8:]
start_time = 0 #defolt 0, s

shot_list = [42121]
for Shotn in shot_list:
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' %Shotn)
    path_res = 'eq_results/%s/' % Shotn
    beta_past = 0

    try:
        os.mkdir(path_res)
    except OSError:
        print('Папка уже существует')

    dia_data = {}
    try:
        p = 0
        dia_file = 'c:/work/equilibrium/dia_data/%d.txt' % Shotn
        if my:
            dia_file = 'c:/work/equilibrium/dia_data/my/%d.txt' % Shotn
        if os.path.exists(dia_file) == False:
            recompliment = int(input('#sht recompliment:'))
            time_step = float(input('time step you need:'))
            get_dia_data.dia_data(Shotn, recompliment, time_step, start_time)
        with open(dia_file, 'r') as file:
            for line in file:
                data = line.split(',')
                if p == 0:
                    for i in data:
                        dia_data[i] = []
                    p += 1
                else:
                    for i, key in enumerate(list(dia_data.keys())):
                        if data[i]:
                            dia_data[key].append(float(data[i]))
                        else:
                            dia_data[key].append(0)
    except:
        print(dia_file)
        print('file not found and file not been created')
        stop
    #print(dia_data)

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

    pdf_file = PdfPages(path_res + tomorow + '_' + str(Shotn) + '_plots.pdf')

    time_list = [i / 1000 for i in dia_data['time']]
    betta_I_sakharov = [round(i, 2) for i in dia_data['betadia\n']]
    betta_I_list = []
    Bt = dia_data['Bt']

    index_start = 0
    # index_start = 23
    index_end = len(time_list)
    # index_end = 24

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
            res_file.write('%14s' % 'alpha1')
            res_file.write('%14s' % 'alpha2')
            res_file.write('\n')

    for ind in range(index_start, index_end):
        print('______________________________________________________________________________________')
        print(ind)
        time = time_list[ind]
        try:
            f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
        except FileNotFoundError:
            print('not found in new version')
            f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
        time, Rc, Rcm, Rav, k, conf = Find_boundary.bound(f, time)
        print('new time: ', time)
        print('configure: ', conf)
        Rav = Rav / 100


        if round(time,4) in time_already_done:
            print('this time alredy done')
            continue
        I_coil = Globus_PET.get_coils(Shotn, time)
        if I_coil_new:
            if I_coil_new[ind] != 0:
                I_coil['Ipl'] = I_coil_new[ind] * 1e3

        for fast_ind in range(9):
            Globus_PET.COIL_upd(Shotn, time, I_coil, Bt[ind], N=fast_ind)
            Globus_PET.DATA_upd(conf, N=fast_ind)

        if beta_past:
            beta_0 = beta_past
        else:
            beta_0 = betta_I_sakharov[ind]

        bound_r = (int(beta_0 * 100) + 10)
        bound_l = (int(beta_0 * 100) - 10)

        li_for_one_time = []
        li_work_one_time = []
        li_work_one_time2 = []
        bI_for_one_time = []
        psi_delta_for_one_time = []
        bounds_delta_for_one_time = []

        psi_delta = 1
        bounds_delta = 1
        attempt = 0
        min_count = 0
        min_delta = 100
        betta_I = beta_0
        beta_min = beta_0 - 0.2
        beta_max = beta_0 + 0.2
        while (abs(psi_delta) > 0.02 or bounds_delta > 0.005) and abs(beta_max - beta_min) > 0.005:
            print('test: ', abs(beta_max - beta_min), abs(beta_max - beta_min) < 0.005)
            '''if attempt % 2:
                betta_I = beta_0 - (attempt%2 + attempt//2) * 0.01
            else:
                betta_I = beta_0 + (attempt % 2 + attempt // 2) * 0.01'''
            betta_I = (beta_max + beta_min) / 2
            print('attempt #%i, beta=%f' %(attempt, betta_I))
            bI_for_one_time.append(betta_I)
            li = 1
            li_f, li_f2, res, psi_delta, bounds_delta, alpha = Globus_PET.find_par_fast(Shotn, time, I_coil, betta_I,
                                                                             pdf_file,
                                                                             psi_dia_sakh=dia_data['psidia'][ind])
            li_for_one_time.append(li_f)
            li_work_one_time.append(alpha)
            li_work_one_time2.append(li_f2)
            psi_delta_for_one_time.append(abs(psi_delta))
            bounds_delta_for_one_time.append(bounds_delta)
            print('res of attempt: ', psi_delta, bounds_delta)
            if psi_delta < 0:
                beta_min = betta_I
            else:
                beta_max = betta_I
            if bounds_delta < min_delta:
                min_delta = bounds_delta
                min_count = 0
            else:
                min_count +=1

            attempt += 1
            if attempt > 30 or min_count > 10:
                break
        if bounds_delta > 0.005:
            if (beta_max + beta_min) / 2 < beta_0:
                beta_max = beta_0 + 0.4
                beta_min = beta_0
            else:
                beta_max = beta_0
                beta_min = beta_0 - 0.4
            while (abs(psi_delta) > 0.02 or bounds_delta > 0.005) and abs(beta_max - beta_min) > 0.005:
                print('test: ', abs(beta_max - beta_min), abs(beta_max - beta_min) < 0.005)
                '''if attempt % 2:
                    betta_I = beta_0 - (attempt%2 + attempt//2) * 0.01
                else:
                    betta_I = beta_0 + (attempt % 2 + attempt // 2) * 0.01'''
                betta_I = (beta_max + beta_min) / 2
                print('attempt #%i, beta=%f' % (attempt, betta_I))
                bI_for_one_time.append(betta_I)
                li = 1
                li_f, li_f2, res, psi_delta, bounds_delta, alpha = Globus_PET.find_par_fast(Shotn, time, I_coil,
                                                                                            betta_I,
                                                                                            pdf_file,
                                                                                            psi_dia_sakh=
                                                                                            dia_data['psidia'][ind])
                li_for_one_time.append(li_f)
                li_work_one_time.append(alpha)
                li_work_one_time2.append(li_f2)
                psi_delta_for_one_time.append(abs(psi_delta))
                bounds_delta_for_one_time.append(bounds_delta)
                print('res of attempt: ', psi_delta, bounds_delta)
                if psi_delta < 0:
                    beta_min = betta_I
                else:
                    beta_max = betta_I
                if bounds_delta < min_delta:
                    min_delta = bounds_delta
                    min_count = 0
                else:
                    min_count += 1

                attempt += 1
                if attempt > 30 or min_count > 10:
                    break

        '''plt.figure()
        plt.plot(psi_delta_for_one_time, 'o-')
        plt.grid()
        plt.ylabel('delta_psi')
        plt.figure()
        plt.plot(li_for_one_time, 'o-')
        plt.grid()
        plt.ylabel('li')
        plt.figure()
        plt.plot(bounds_delta_for_one_time, 'o-')
        plt.grid()
        plt.ylabel('delta_bounds')
        plt.figure()
        plt.plot(bI_for_one_time, 'o-')
        plt.grid()
        plt.ylabel('beta_I')
        plt.show()'''

        ind_min = bounds_delta_for_one_time.index(min(bounds_delta_for_one_time))
        print('before check')
        print(bounds_delta_for_one_time[ind_min], psi_delta_for_one_time[ind_min])
        mera_list1 = []
        mera_list2 = []
        if psi_delta_for_one_time[ind_min] > 0.04:
            for i in range(len(bounds_delta_for_one_time)):
                mera_list1.append((psi_delta_for_one_time[i]/10 + bounds_delta_for_one_time[i]*10) / 2)
                mera_list2.append((psi_delta_for_one_time[i] ** 2 + bounds_delta_for_one_time[i] ** 2) ** 0.5)
            print(min(mera_list1), mera_list1.index(min(mera_list1)),
                  bounds_delta_for_one_time[mera_list1.index(min(mera_list1))],
                  psi_delta_for_one_time[mera_list1.index(min(mera_list1))])
            print(min(mera_list2), mera_list2.index(min(mera_list2)),
                  bounds_delta_for_one_time[mera_list2.index(min(mera_list2))],
                  psi_delta_for_one_time[mera_list2.index(min(mera_list2))])
            ind_min = mera_list1.index(min(mera_list1))
        print('after check')
        print(bounds_delta_for_one_time[ind_min], psi_delta_for_one_time[ind_min])

        betta_I = bI_for_one_time[ind_min]
        li_acc2 = li_work_one_time[ind_min]

        b_I, li_3, bp, W_all, We, V, S, P_axis, li_code, Ftor_pl, Wi, bounds_delta_res, ne_av, strike_point, q95 = Globus_PET.find_bound(
            Shotn, time, I_coil, betta_I, li_acc2, 0, 1, pdf=pdf_file,
            inside=True, Wi=True, share=True)

        li_code = li_work_one_time2[ind_min]
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
        W_list.append(betta_I * (mu0 * I_coil['Ipl'] * 1000 * I_coil['Ipl'] * 1000 * Rcm) / 4)
        with open(path_res + 'test_' + tomorow + 'res_' + str(Shotn) + '.txt', 'a') as res_file:
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
            res_file.write('%8.4f' % float(li_code))
            res_file.write('%14.4f' % Wi)
            res_file.write('%8.4f' % float(bounds_delta_for_one_time[ind_min]))
            res_file.write('%8.4f' % float(bounds_delta_res))
            res_file.write('%8.2f' % I_coil['Ipl'])
            res_file.write('%8.4f' % Bt[ind])
            res_file.write('%14.4e' % ne_av)
            res_file.write('%14.4f' % li_acc2[0])
            res_file.write('%14.4f' % li_acc2[1])
            res_file.write('\n')


        beta_past = float(betta_I)
    pdf_file.close()
    print('time left: ', - start_time + t.time())

