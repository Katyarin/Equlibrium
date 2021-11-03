import niifaRipper
import matplotlib.pyplot as plt
import ripper
import numpy as np
import time as timer
import numpy as np
from scipy import interpolate as interp
import json

shotn = [i for i in range(37700, 40600, 29)]


def linear_approx(time_arr, sig_arr, time_ind, x):
    x1 = time_arr[time_ind[0]]
    x2 = time_arr[time_ind[1]]
    y1 = sig_arr[time_ind[0]]
    y2 = sig_arr[time_ind[1]]
    return y2 + (y1 - y2) / (x1 - x2) * (x - x2)


def find_sig(i, shot):
    signals = niifaRipper.extract_niifa('//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/magn_data/', str(shot) + '_' + str(i))
    if abs(signals['Ipl'][int(len(signals['Ipl']) / 2)]) < 50:
        find_sig(i+1, shot)
    else:
        print(i, 'success!')


def get_data(shotn, data_name):
    raw_data, link = ripper.extract('//172.16.12.127/data', shotn, data_name)
    data = {}
    for key in dict.keys(raw_data):
        data[key] = {}
        data[key]['name'] = raw_data[key]['name']
        data[key]['time'], data[key]['data'] = ripper.x_y(raw_data[key])

    data_clean = {}
    data_clean['Shotn'] = shotn
    #data_clean['time'] = data[list(dict.keys(data))[0]]['time']

    for key in dict.keys(data):
        data_clean[data[key]['name']] = {}
        data_clean[data[key]['name']]['data'] = data[key]['data']
        data_clean[data[key]['name']]['time'] = data[key]['time']

    comb_change_name = ['Shotn', 'Ipl', 'Ics', 'Ihfc', 'Ivfc', 'Icc', 'Ipf1', 'Ipf3']
    keys = list(data_clean.keys())
    print(data_clean.keys())
    beauty_data = {}

    for key in range(len(data_clean.keys())):
        beauty_data[comb_change_name[key]] = data_clean[keys[key]]

    return beauty_data


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


i_need = ['Ipf1', 'Ipf3', 'Ics', 'Icc', 'Ihfc', 'Ivfc', 'Ipl']
all_diff = {}
all_diff['globus'] = {}
all_diff['niifa'] = {}
for I in i_need:
    all_diff['globus'][I] = []
    all_diff['niifa'][I] = []
all_diff['time'] = []
all_diff['shots'] = []

start_time = timer.time()
for shot in shotn:
    try:
        signals = niifaRipper.extract_niifa('//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/magn_data/', str(shot))
    except FileNotFoundError:
        continue
    if abs(signals['Ipl'][int(len(signals['Ipl']) / 2)]) < 50:
        try:
            find_sig(0, shot)
        except FileNotFoundError:
            continue
    print(shot)
    name_Ip_in_comb = ['Ip внутр', 'Ipf1', 'Ipf3', 'Ihfc', 'Ivfc', 'Icc', 'Ics']
    combiscope_data = get_data(shot, name_Ip_in_comb)

    all_diff['shots'].append(shot)
    #all_diff['time'].append(signals['t_ms'])

    signals_smooth = {}
    niifa = {}
    globus = {}
    combiscope_data_smooth = {}
    for I in i_need:
        niifa[I] = []
        globus[I] = []
        signals_smooth[I] = smooth(signals[I], 16)
        if I != 'Ipl':
            combiscope_data_smooth[I] = smooth([i / 1000 for i in combiscope_data[I]['data']], 100)
        else:
            combiscope_data_smooth[I] = smooth([i / 1000 for i in combiscope_data[I]['data']], 1000)

    if signals['t_ms'][-1] < combiscope_data[I]['time'][-1] * 1000:
        last_ms = int(signals['t_ms'][-1])
    else:
        last_ms = int(combiscope_data[I]['time'][-1] * 1000)
    print(last_ms)
    all_diff['time'].append([i for i in range(1, last_ms)])

    for I in i_need:
        sig_interp = interp.interp1d(signals['t_ms'], signals_smooth[I], kind='linear')
        comb_interp = interp.interp1d([i * 1000 for i in combiscope_data[I]['time']], combiscope_data_smooth[I],
                                      kind='linear')
        for ms in range(2, last_ms):
            niifa[I].append(float(sig_interp(ms)))
            try:
                globus[I].append(float(comb_interp(ms)))
            except ValueError:
                print(ms, len(combiscope_data[I]['time']), combiscope_data[I]['time'][0], combiscope_data[I]['time'][-1])
                stop
        all_diff['globus'][I].append(globus[I])
        all_diff['niifa'][I].append(niifa[I])
print(timer.time() - start_time)
with open('shot_compare_100.json', 'w') as file:
    json.dump(all_diff, file)






