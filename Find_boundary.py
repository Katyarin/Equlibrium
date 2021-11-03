import json
import matplotlib.pyplot as plt

f = 'mcc_37893.json'

def bound(file_name, t, plotting=False, ax=0):
    with open(file_name, 'r') as file:
        data_mcc = json.load(file)
        time = data_mcc['time']['variable']
        rb = data_mcc['boundary']['rbdy']['variable']
        zb = data_mcc['boundary']['zbdy']['variable']
        data_mcc = []
    index_time = 0
    for i in range(len(time)):
        if time[i] > t - 0.0005  and time[i] < t + 0.0005:
            print( t, t - 0.0005, t + 0.0005)
            print('MCC time: ', time[i])
            index_time = i
    if index_time == 0:
        print('Index time not found')
        stop
        return 0
    boundary = {'time': time[index_time], 'r': rb[index_time][1:], 'z': zb[index_time][1:]}
    Rc = (max(rb[index_time][1:]) + min(rb[index_time][1:])) / 200
    Rcm = sum(rb[index_time][1:]) / len(rb[index_time][1:]) / 100
    print(max(rb[index_time][1:]), min(rb[index_time][1:]))
    with open('c:/work/equilibrium/Globus_PET/MCC.json', 'w') as file2:
        json.dump(boundary, file2)
    if plotting == True:
        print('we_here')
        ax.plot(boundary['r'], boundary['z'])
    return time[index_time], Rc, Rcm


