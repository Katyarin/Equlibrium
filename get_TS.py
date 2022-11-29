import requests
import json
import matplotlib.pyplot as plt
URL = 'https://172.16.12.130:443/api'
q = 1.602176634e-19
PATH = {22: '//172.16.12.127/Pub/!!!TS_RESULTS/2022/', 21:  '//172.16.12.127/Pub/!!!TS_RESULTS/2021/',
        20: '//172.16.12.127/Pub/!!!TS_RESULTS/2020/'}

'''shotn = 41615
time = 180'''

def get_data(shotn, time):
    ne_raw = []
    Te_raw = []
    for year in PATH.keys():
        try:
            print(PATH[year] + str(shotn) + '/' +str(shotn) + '_n(R).csv')
            with open(PATH[year] + str(shotn) + '/' +str(shotn) + '_n(R).csv', 'r') as nfile:
                for line in nfile:
                    ne_raw.append(line.split(','))
            with open(PATH[year] + str(shotn) + '/' + str(shotn) + '_T(R).csv', 'r') as nfile:
                for line in nfile:
                    Te_raw.append(line.split(', '))
            error = None
            break
        except:
            print('not find in ' + str(year))
            continue

    if ne_raw:
        for i in range(1, len(ne_raw[0]), 2):
            if time - 2 < float(ne_raw[0][i]) < time + 2:
                print('TS_time = ', float(ne_raw[0][i]))
                time_ind = i
        R = [float(ne_raw[i][0]) for i in range(2, len(ne_raw)) if Te_raw[i][time_ind] != '--']
        Te = [float(Te_raw[i][time_ind]) for i in range(2, len(Te_raw)) if Te_raw[i][time_ind] != '--']
        Te_err = [float(Te_raw[i][time_ind+1]) for i in range(2, len(Te_raw)) if Te_raw[i][time_ind] != '--']
        ne = [float(ne_raw[i][time_ind]) for i in range(2, len(ne_raw)) if Te_raw[i][time_ind] != '--']
        ne_err = [float(ne_raw[i][time_ind+1]) for i in range(2, len(ne_raw)) if Te_raw[i][time_ind] != '--']
        print('All find')
    else:
        error = 'Mb not really ne'
        print('not find in files')
        response = requests.post(url=URL, verify=False, json={
            'subsystem': 'db',
            'reqtype': 'get_shot_verified',
            'shotn': int(shotn)
        })
        try:
            data = response.json()
            with open('dump.json', 'w') as file:
                json.dump(data, file)
        except:
            print('Not a json?')
            response = requests.post(url=URL, verify=False, json={
                'subsystem': 'db',
                'reqtype': 'get_shot',
                'shotn': int(shotn)
            })
            data = response.json()
            with open('dump.json', 'w') as file:
                json.dump(data, file)

        for event in data['data']['events']:
            if event['error'] != None:
                continue
            if time - 2 < event['timestamp'] < time + 2:
                print('TS_time = ', event['timestamp'])
                try:
                    R = [i['R'] for i in data['data']['config']['poly']]
                except KeyError:
                    try:
                        R = [i['R'] for i in data['data']['polys']]
                    except KeyError:
                        R = [i['R'] for i in data['data']['config']['fibers']]
                        print(R)
                Te = [event['T_e'][poly]['T'] for poly in range(10)]
                Te_err = [event['T_e'][poly]['Terr'] for poly in range(10)]
                ne = [event['T_e'][poly]['n'] for poly in range(10)]
                ne_err = [event['T_e'][poly]['n_err'] for poly in range(10)]
                print('All find')

    P = [Te[i]*ne[i]*q for i in range(len(R))]
    P_err = [((Te_err[i]/Te[i])**2 + (ne_err[i]/ne[i])**2)**0.5 * P[i] for i in range(len(R))]
    '''print(Te)
    print(Te_err)
    print(ne)
    print(ne_err)
    print(P)
    print(P_err)'''

    return R, {'Te': Te, 'err': Te_err}, {'ne': ne, 'err': ne_err}, {'P': P, 'err': P_err}, error