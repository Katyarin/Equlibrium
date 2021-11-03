import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

def li_approx(alf1, alf2):
    a = -0.049 + 1.34 * alf2
    b1 = 0.41 - 0.944 * alf2
    b2 = -0.1345 + 0.26 * alf2
    return a + b1 * alf1 + b2 * alf1 * alf1

def find_li(li):
    alf1 = (li - 0.2 + 0.6) / 2
    a1 = -0.049
    a2 = 1.34
    b11 = 0.41
    b12 = - 0.944
    b21 = -0.1345
    b22 = 0.26
    alf2 = (li - a1 - b11 * alf1 - b21 * alf1 * alf1) / (a2 + b12 * alf1 + b22 * alf1 * alf1)
    return alf1, alf2

Shotn = 38800
time_list = []
W_list = []
betta_p_list = []
li_find_list = []
with open('W' + str(Shotn) + '.txt', 'r') as wfile:
    for line in wfile:
        time_list.append(float(line.split()[0]))
        W_list.append(float(line.split()[1]))
'''W_list=[900]
time_list = [0.2]'''

#for ind in range(len(time_list[10:11])):
for ind in [17, 18]:
    print('______________________________________________________________________________________')
    print(ind)
    time = time_list[ind]
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/new_mcc/mcc_%d.json' % Shotn
    time, Rc = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    Globus_PET.COIL_upd(Shotn, time, I_coil)

    betta_pol = 8/3 * W_list[ind] / (mu0*I_coil['Ipl']*1000*I_coil['Ipl']*1000*Rc)
    print('betta_pol: ', betta_pol)
    betta_p_list.append(betta_pol)

    betta_I = betta_pol / 1.2
    print('betta_I: ', betta_I)
    li = 1
    #alf11, alf22 = 1.5, 1
    alf11, alf22 = find_li(1.3)
    print(alf11, alf22)
    print('li_by_approx: ', li_approx(alf11, alf22))
    #li_find, res = Globus_PET.find_par('li', Shotn, time, I_coil, betta_I, li, [80, 140, 5], show2=True)
    Globus_PET.DURS_upd(Shotn, time, I_coil['Ipl'], betta_I, alf11, alf22)
    try:
        process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(timeout=5)
        # print('out: ', stdout.decode('utf-8'))
        result = Globus_PET.EVOL_res()
        if result == -1:
            print('NOT COUNT')
            continue
        dif_x, dif_y = Globus_PET.compare_bound(0, show=True)
        print('real li: ', result['li'])
        print('real bp: ', result['betpol'])
        Rc2 = result['Rc']
        print(Rc, I_coil['Ipl'])
        print(Rc2)
        print('~bpolin = ', result['li'] * Rc2 * (I_coil['Ipl'])**2 / (4 * pi))
        plt.show()
    except subprocess.TimeoutExpired:
        print('time over')
        subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))

