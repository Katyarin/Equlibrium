import ripper
import numpy as np
import niifaRipper
PATH = 'c:/work/equilibrium/Globus_PET/'
COIL = []


'''def get_data(shotn, data_name):
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

    comb_change_name = ['Shotn', 'Ip', 'CS', 'HFC', 'VFC', 'CC', 'PF1', 'PF3']
    keys = list(data_clean.keys())
    print(data_clean.keys())
    beauty_data = {}

    for key in range(len(data_clean.keys())):
        beauty_data[comb_change_name[key]] = data_clean[keys[key]]

    return beauty_data'''

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

with open(PATH + 'Start_files/COIL.DAT', 'r') as file:
    for line in file:
        COIL.append(line)

Shotn = int(COIL[0][39:44])
time = float(COIL[0][47:51]) #sec
Ip = float(COIL[0][57:62]) #MA

n_str_coil = [3, 8, 13, 18, 27, 32, 45]
coil_name = ['Ipf1', 'Ipf2', 'Ipf3', 'Ihfc', 'Ivfc', 'Icc', 'Ics']
I_coil = {}

for i in range(len(n_str_coil)):
    if n_str_coil[i] < 40:
        I_coil[coil_name[i]] = float(COIL[n_str_coil[i]][63:78]) * 1e6
    else:
        I_coil[coil_name[i]] = float(COIL[n_str_coil[i]][55:68]) * 1e6

DURS = []
with open(PATH + 'Start_files/DURS.DAT', 'r') as file2:
    for line in file2:
        DURS.append(line)

B_pol = float(DURS[1][1:6]) #MA
Ip = float(DURS[2][1:7]) #MA

I_coil['Ipl'] = Ip * 1e6

print(I_coil)

#name_Ip_in_comb = ['Ip внутр', 'Ipf1', 'Ipf3', 'Ihfc', 'Ivfc', 'Icc', 'Ics']

#combiscope_data = get_data(Shotn, name_Ip_in_comb)

#print(combiscope_data.keys())
signals = niifaRipper.extract_niifa('//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/magn_data/', Shotn)

i_need = ['Ipf1', 'Ipf3', 'Ics', 'Ipl', 'Icc', 'Ihfc', 'Ivfc', 'Ipf2']

signals['Ipf2'] = [0.001] * len(signals['t_ms'])

I_coil_comb = {}
k = 0
for I in i_need:
    k = 0
    for t_mom in range(len(signals[I])):
        if signals['t_ms'][t_mom] / 1000 > time - 0.00004 and signals['t_ms'][t_mom] / 1000 < time + 0.00004:
            index_time2 = t_mom
            k+=1
            continue
    if k > 1:
        print('many times')
        stop
    if I == 'Ihfc':
        I_coil_comb[I] = smooth(signals[I], 95)[index_time2] * 1e3 #kA
    else:
        I_coil_comb[I] = signals[I][index_time2] * 1e3 #kA
Difference = {}

'''combiscope_data['PF2'] = {}
combiscope_data['PF2']['time'] = [time]
combiscope_data['PF2']['data'] = [1]'''

for I in I_coil.keys():
    print(I)
    print(I_coil[I], I_coil_comb[I])
    Difference[I] = I_coil[I] / I_coil_comb[I]

print(Difference)
