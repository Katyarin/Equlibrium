import numpy as np
import Find_boundary
import matplotlib.pyplot as plt

Shotn = 38800

time_list = []
W_list = []
with open('W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))

f1 = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/new_mcc/mcc_%d.json' % Shotn
f2 = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_data/mcc_%d.json' % Shotn
time = 0.145

for time in time_list:
    fig = plt.figure(figsize=(5, 8))
    ax = fig.add_subplot(111)
    ax.grid()
    time1, Rc1 = Find_boundary.bound(f1, time, plotting=True, ax=ax)
    time2, Rc2 = Find_boundary.bound(f2, time, plotting=True, ax=ax)
    print(time1, time2)
    print(Rc1, Rc2)
    plt.show()