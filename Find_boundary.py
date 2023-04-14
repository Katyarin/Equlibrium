import json
import matplotlib.pyplot as plt

f = 'mcc_37893.json'

def bound(file_name, t, plotting=False, ax=0):
    lim_min = 0.125
    lim_max = 0.61
    with open(file_name, 'r') as file:
        data_mcc = json.load(file)
        time = data_mcc['time']['variable']
        rb = data_mcc['boundary']['rbdy']['variable']
        zb = data_mcc['boundary']['zbdy']['variable']
        r_c = data_mcc['current_coils']['r']['variable']
        I_c = data_mcc['current_coils']['I']['variable']
        Ip = data_mcc['Ipl']['variable']
        data_mcc = []
    index_time = 0
    for i in range(len(time)):
        if time[i] > t - 0.0005  and time[i] < t + 0.0005:
            #print( t, t - 0.0005, t + 0.0005)
            #print('MCC time: ', time[i])
            index_time = i
    if index_time == 0:
        print('Index time %f not found' %t)
        stop
        return 0
    boundary = {'time': time[index_time], 'r': rb[index_time][1:], 'z': zb[index_time][1:]}
    Rav = sum([r_c[index_time][i] * I_c[index_time][i] for i in range(len(r_c[index_time]))]) / Ip[index_time]
    k = (max(boundary['z']) - min(boundary['z'])) / (max(boundary['r']) - min(boundary['r']))
    Rc = (max(rb[index_time][1:]) + min(rb[index_time][1:])) / 200
    Rcm = sum(rb[index_time][1:]) / len(rb[index_time][1:]) / 100
    if abs(max(boundary['r']) / 100 - lim_max) < 0.001 or abs(min(boundary['r']) / 100 - lim_min) < 0.001:
        print('low field side', max(boundary['r']) / 100 - lim_max)
        print('high field side', min(boundary['r']) / 100 - lim_min)
        configure = 'lim'
    else:
        print('low field side', max(boundary['r']) / 100 - lim_max)
        print('high field side', min(boundary['r']) / 100 - lim_min)
        configure = 'div'
    #print(max(rb[index_time][1:]), min(rb[index_time][1:]))
    with open('c:/work/equilibrium/Globus_PET/MCC.json', 'w') as file2:
        json.dump(boundary, file2)
    if plotting == True:
        print('we_here')
        ax.plot(boundary['r'], boundary['z'])
    return time[index_time], Rc, Rcm, Rav, k, configure


