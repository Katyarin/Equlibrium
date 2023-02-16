import shtRipper
import Find_boundary
import matplotlib.pyplot as plt
import numpy as np
import json
pi = 3.14159265359

Globus3 = False

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_sht_data(Shotn, data_name):
    filename = '//172.16.12.127/Data/sht%i.sht' %Shotn
    if Globus3:
        filename = '//172.16.12.28/Data/sht%i.sht' % Shotn
    res = shtRipper.ripper.read(filename, data_name)
    smooth_res = {}
    print(res.keys())
    for key in res.keys():
        if key != 'nl 42 cm (1.5мм) 64pi':
            baseline = sum(res[key]['y'][:1000]) / len(res[key]['y'][:1000])
            res[key]['y'] = [i - baseline for i in res[key]['y']]
        smooth_res[key] = {}
        if key == 'Ip внутр.(Пр2ВК) (инт.18)':
            smooth_res[key] = {}
            smooth_res[key]['data'] = list(smooth(res[key]['y'][5:-1:10], 31))
            smooth_res[key]['time'] = res[key]['x'][5:-1:10]
        else:
            smooth_res[key]['data'] = list(smooth(res[key]['y'], 95))
            smooth_res[key]['time'] = res[key]['x']
    return smooth_res

def dia_data(shot, recoupment, delta_time, start_time=0, end_time=0):

    data_name_need = ['Ip внутр.(Пр2ВК) (инт.18)', 'Itf (2TF)(инт.23)', 'Диамагнитный сигнал (новый инт.)', 'Ics (4CS) (инт.22)',
                  'Up (внутреннее 175 петля)', 'Программа тока Ip']


    data = get_sht_data(shot, data_name_need)
    recoupment_data = get_sht_data(recoupment, data_name_need)

    a = 1
    t_start = 0.14
    t_end = 0.24
    if 'Программа тока Ip' in list(data.keys()):
        programm_Ip_diff = [(data['Программа тока Ip']['data'][i+1] - data['Программа тока Ip']['data'][i]) / (data['Программа тока Ip']['time'][i+1] - data['Программа тока Ip']['time'][i]) for i in range(len(data['Программа тока Ip']['time']) - 1)]
        for i, t in enumerate(data['Программа тока Ip']['time'][:-1]):
            if 0.24 > t > 0.13:
                if programm_Ip_diff[i] < 20 * a:
                    if a > 0:
                        t_start = round(t, 2)
                        a = -a
                    else:
                        t_end = round(t, 2)
                        break
        if start_time:
            t_start = start_time
        if end_time:
            t_end = end_time
        plt.figure()
        plt.plot(data['Программа тока Ip']['time'], data['Программа тока Ip']['data'])
        plt.axvline(t_start, color='r')
        plt.axvline(t_end, color='r')
        plt.grid()

        plt.figure()
        plt.plot(data['Программа тока Ip']['time'][:-1], programm_Ip_diff, 'o-')
        plt.axvline(t_start, color='r')
        plt.axvline(t_end, color='r')
    print(t_start, t_end)
    time_need = [i/10 for i in range(int((t_start*10000)+ delta_time*10), int(t_end*10000), int(delta_time*10))]
    print(time_need)
    for i in [1,2]:
        plt.figure()
        plt.title(data_name_need[i])
        plt.plot(data[data_name_need[i]]['time'], data[data_name_need[i]]['data'], label=shot)
        plt.plot(recoupment_data[data_name_need[i]]['time'], recoupment_data[data_name_need[i]]['data'], label='recoupment #' + str(recoupment))
        plt.grid()
        plt.axvline(t_start, color='r')
        plt.axvline(t_end, color='r')
        plt.legend()

    plt.figure()
    plt.plot(data['Ip внутр.(Пр2ВК) (инт.18)']['time'], data['Ip внутр.(Пр2ВК) (инт.18)']['data'])
    plt.axvline(t_start, color='r')
    plt.axvline(t_end, color='r')
    plt.show()


    dia_sig1 = [data['Диамагнитный сигнал (новый инт.)']['data'][i] + data['Ics (4CS) (инт.22)']['data'][i] * 8e-6 for i in range(len(data['Ics (4CS) (инт.22)']['data']))]
    dia_sig2 = [recoupment_data['Диамагнитный сигнал (новый инт.)']['data'][i] + recoupment_data['Ics (4CS) (инт.22)']['data'][i] * 8e-6 for i in range(len(data['Ics (4CS) (инт.22)']['data']))]

    diamagnetic_sig = {'time': data['Диамагнитный сигнал (новый инт.)']['time'],
                       'data': [(dia_sig1[i] - dia_sig2[i]) * 2.915 for i in range(len(dia_sig1))]}

    '''plasma current'''
    Ip_all = [data['Ip внутр.(Пр2ВК) (инт.18)']['data'][i] - data['Up (внутреннее 175 петля)']['data'][i] * 370 - recoupment_data['Ip внутр.(Пр2ВК) (инт.18)']['data'][i] for i in range(len(data['Ip внутр.(Пр2ВК) (инт.18)']['data']))]

    '''Bt'''
    Rc_list =[]
    Rcm_list = []
    Rav_list = []
    k_list = []
    for t in time_need:
        try:
            f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % shot
        except FileNotFoundError:
            print('not found in new version')
        # f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
        time, Rc, Rcm, Rav, k = Find_boundary.bound(f, t/1000)
        Rc_list.append(Rc)
        Rcm_list.append(Rcm)
        Rav_list.append(Rav/100)
        k_list.append(k)
    with open('c:/work/equilibrium/dia_data/my/%s.txt' %shot, 'w') as file:
        file.write('time,')
        file.write('psidia,')
        file.write('Ip,')
        file.write('Bt,')
        file.write('betadia')
        file.write('\n')
        for i, t in enumerate(diamagnetic_sig['time']):
            for j, p in enumerate(time_need):
                if p + 0.006 > t*1000 > p - 0.001 :
                    Bt = 0.2* 16 * data['Itf (2TF)(инт.23)']['data'][i] /1e6/Rav_list[j]
                    beta_dia = 1- (k_list[j] *k_list[j] +1) / (2*k_list[j]) * Bt *diamagnetic_sig['data'][i] / (20*pi*Ip_all[i]/1e6*Ip_all[i]/1e6)
                    print(Ip_all[i] / 1e6, Bt, beta_dia)
                    #print(p, diamagnetic_sig['data'][i], Ip_all[i]/1e6, Bt, beta_dia)
                    file.write('%4.1f,' %p)
                    file.write('%5.4f,' % diamagnetic_sig['data'][i])
                    file.write('%5.4f,' % (Ip_all[i]/1e6))
                    file.write('%5.4f,' % Bt)
                    file.write('%5.4f' % beta_dia)
                    file.write('\n')
        print('dia_file saved')


shotn = 42770
rec = 42776
dia_data(shotn, rec, 2, end_time=0)