import requests
import json
import matplotlib.pyplot as plt
URL = 'https://172.16.12.130:443/api'
q = 1.602176634e-19

'''shotn = 41615
time = 180'''

def get_data(shotn, time):
    response = requests.post(url=URL, verify=False, json={
        'subsystem': 'db',
        'reqtype': 'get_shot',
        'shotn': int(shotn)
    })
    try:
        data = response.json()
        with open('dump.json', 'w') as file:
            json.dump(data, file)
    except:
        print('Not a json?')

    for event in data['data']['events']:
        if event['error'] != None:
            continue
        if time - 2 < event['timestamp'] < time + 2:
            R = [i['R'] for i in data['data']['config']['poly']]
            Te = [event['T_e'][poly]['T'] for poly in range(10)]
            Te_err = [event['T_e'][poly]['Terr'] for poly in range(10)]
            ne = [event['T_e'][poly]['n'] for poly in range(10)]
            ne_err = [event['T_e'][poly]['n_err'] for poly in range(10)]
            print('All find')

    P = [Te[i]*ne[i]*q for i in range(len(R))]
    P_err = [(Te_err[i]**2 + ne_err[i]**2)**0.5 * q for i in range(len(R))]

    return R, {'Te': Te, 'err': Te_err}, {'ne': ne, 'err': ne_err}, {'P': P, 'err': P_err}