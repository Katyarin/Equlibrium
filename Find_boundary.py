import json

f = 'mcc_37893.json'

def bound(file_name, t):
    with open(file_name, 'r') as file:
        data_mcc = json.load(file)
        time = data_mcc['time']['variable']
        rb = data_mcc['boundary']['rbdy']['variable']
        zb = data_mcc['boundary']['zbdy']['variable']
        data_mcc = []
    index_time = 0
    for i in range(len(time)):
        if time[i] > t - 0.0001  and time[i] < t + 0.0001:
            print('MCC time: ', time[i])
            index_time = i
    if index_time == 0:
        print('Index time not found')
        stop
        return 0
    boundary = {'time': time[index_time], 'r': rb[index_time], 'z': zb[index_time]}
    with open('c:/work/equilibrium/Globus_PET/MCC.json', 'w') as file2:
        json.dump(boundary, file2)
    return time[index_time]


