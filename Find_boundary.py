import json
import matplotlib.pyplot as plt

f = 'mcc_37893.json'

def bound(file_name, t, plotting=False, ax=0, beta_li=0, Psi_av=0, tr=0, xDot=0):
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
        if beta_li:
            beta_li_list = data_mcc['shafr_int_method']['beta_eq+li/2']['variable']
            Bzav = data_mcc['Bzav']['variable']
        if Psi_av:
            Psav = data_mcc['Psav']['variable']
    index_time = 0
    for i in range(len(time)):
        if time[i] > t - 0.0005  and time[i] < t + 0.0005:
            #print( t, t - 0.0005, t + 0.0005)
            #print('MCC time: ', time[i])
            index_time = i
    if index_time == 0:
        print('Index time %f not found' %t)
        return 0
    boundary = {'time': time[index_time], 'r': rb[index_time][1:], 'z': zb[index_time][1:]}
    Rav = sum([r_c[index_time][i] * I_c[index_time][i] for i in range(len(r_c[index_time]))]) / Ip[index_time]
    k = (max(boundary['z']) - min(boundary['z'])) / (max(boundary['r']) - min(boundary['r']))
    Rc = (max(rb[index_time][1:]) + min(rb[index_time][1:])) / 200
    Rcm = sum(rb[index_time][1:]) / len(rb[index_time][1:]) / 100
    a = (max(rb[index_time][1:]) - min(rb[index_time][1:])) / 200
    if beta_li:
        beta_li = beta_li_list[index_time]
        Bv = Bzav[index_time]
    if Psi_av:
        Psav_local = Psav[index_time]
    if tr:
        z_up = max(zb[index_time][1:])
        z_low = min(zb[index_time][1:])
        r_up = rb[index_time][zb[index_time].index(z_up)] / 100
        r_low = rb[index_time][zb[index_time].index(z_low)] /100
        tr_up = (Rc - r_low) / a
        tr_down = (Rc - r_up) / a
    if abs(max(boundary['r']) / 100 - lim_max) < 0.001 or abs(min(boundary['r']) / 100 - lim_min) < 0.001:
        #print('low field side', max(boundary['r']) / 100 - lim_max)
        #print('high field side', min(boundary['r']) / 100 - lim_min)
        configure = 'lim'
    else:
        #print('low field side', max(boundary['r']) / 100 - lim_max)
        #print('high field side', min(boundary['r']) / 100 - lim_min)
        configure = 'div'
    print(configure)
    if xDot:
        xCords = (data_mcc['Rx']['variable'][index_time]/100, data_mcc['Zx']['variable'][index_time]/100)
    #print(max(rb[index_time][1:]), min(rb[index_time][1:]))
    with open('c:/work/equilibrium/Globus_PET/MCC.json', 'w') as file2:
        json.dump(boundary, file2)
    with open('results/%s_%.3f.txt' %(file_name[-10:-5], t), 'w') as file2:
        file2.write('%10s' %'r')
        file2.write('%10s' %'z')
        file2.write('\n')
        for i in range(len(boundary['r'])):
            file2.write('%10.3f' %boundary['r'][i])
            file2.write('%10.3f' % boundary['z'][i])
            file2.write('\n')
    if plotting == True:
        print('we_here')
        ax.plot([i/100 for i in boundary['r']], [i/100 for i in boundary['z']])
    if beta_li:
        return time[index_time], Rc, Rcm, Rav, k, configure, beta_li, a, Bv
    if Psi_av:
        return time[index_time], Psav_local
    if tr:
        return time[index_time], Rc, Rcm, Rav, k, configure, tr_up, tr_down, Ip[index_time]
    if xDot:
        return time[index_time], Rc, Rcm, Rav, k, configure, xCords
    return time[index_time], Rc, Rcm, Rav, k, configure

'''Shotn = 43115

PATH = 'c:/work/equilibrium/Globus_PET/'
with open(PATH + 'BLANFW.DAT', 'r') as file3:
    blanfw = []
    for line in file3:
        blanfw.append(line.split())

Rvv = float(blanfw[0][0])
nvv = int(blanfw[1][0])

with open(PATH + 'LIMPNT.dat', 'r') as file4:
    limpnt = []
    for line in file4:
        limpnt.append(line.split())

nlim = int(limpnt[0][0])

f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
for t in range(120, 220, 10):
    fig = plt.figure(figsize=(5, 8))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, 0.8)
    bound(f, t/1000, plotting=True, ax=ax)
    ax.set_title(t)
    ax.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
    ax.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')

    ax.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
            [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')
plt.show()'''



