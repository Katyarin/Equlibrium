import copy
import random
import subprocess

from matplotlib import pyplot as plt

import Find_boundary
import Globus_PET
import json
PATH = 'c:/work/equilibrium/Globus_PET/'

def linInter(xlist:list, ylist:list, xNew):
    for i in range(len(xlist)-1):
        if xlist[i] <= xNew <= xlist[i+1] or xlist[i] >= xNew >= xlist[i+1]:
            return ylist[i] + (ylist[i+1] - ylist[i]) * (xNew - xlist[i]) / (xlist[i+1] - xlist[i])
    return None

def findZ(xlist:list, ylist:list, x1, y1):
    difList = [((xlist[i] - x1)**2 + (ylist[i] - y1)**2)**0.5 for i in range(len(xlist))]
    return ylist[difList.index(min(difList))]

def chi2_bound(boundMCC, boundPET):
    upBoundMCC = {}
    upBoundMCC['zb'] = [el for i, el in enumerate(boundMCC['zb']) if boundMCC['zb'][i] >= 0]
    upBoundMCC['rb'] = [el for i, el in enumerate(boundMCC['rb']) if boundMCC['zb'][i] >= 0]
    downBoundMCC = {}
    downBoundMCC['zb'] = [el for i, el in enumerate(boundMCC['zb']) if boundMCC['zb'][i] < 0]
    downBoundMCC['rb'] = [el for i, el in enumerate(boundMCC['rb']) if boundMCC['zb'][i] < 0]
    chi2 = 0
    for i, r in enumerate(boundPET['rb']):
        zMCC = findZ(boundMCC['rb'], boundMCC['zb'], r, boundPET['zb'][i])
        chi2 += (zMCC - boundPET['zb'][i]) * (zMCC - boundPET['zb'][i])
    return chi2

def MCCbound(Shotn, time):
    try:
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
    except FileNotFoundError:
        print('not found in new version')
        f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    with open(f, 'r') as file:
        data_mcc = json.load(file)
    timeMCC = data_mcc['time']['variable']
    timeIndex = 0
    for i in range(len(timeMCC)-1):
        if timeMCC[i+1] > time*1e-3 >= timeMCC[i]:
            if abs(timeMCC[i+1] - time*1e-3) > abs(time*1e-3 - timeMCC[i]):
                timeIndex = i
            else:
                timeIndex = i + 1
    MCC = {'rb': [i/100 for i in data_mcc['boundary']['rbdy']['variable'][timeIndex]],
           'zb': [i/100 for i in data_mcc['boundary']['zbdy']['variable'][timeIndex]]}

    return MCC

def runGlobusPET(Shotn, time, Ipl, betta_I, alf11, alf22, N=0):
    Globus_PET.DURS_upd(Shotn, time, Ipl, betta_I, alf11, alf22, N)
    try:
        process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdout, stderr = process.communicate(timeout=5)
        result = Globus_PET.EVOL_res()
        if result == -1:
            print('NOT COUNT')
            return -1
    except subprocess.TimeoutExpired:
        print('time over')
        subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))
        return -1
    return result


def AvMinMethod(a, b):
    par1 = (3 * a + b) / 4
    par2 = (3 * b + a) / 4
    return par1, par2


def getPetBound(Shotn, time, N=0):
    if N:
        PATH_loc = PATH + 'PET' + str(N) + '/'
    else:
        PATH_loc = PATH

    rgr, zgr, ugr, up, rx0, zx0, rm, zm, norm_ugr, ipr, js0 = Globus_PET.out_data_read(Shotn, time, PATH_loc, share=0)
    cs = plt.contour(rgr, zgr, ugr, levels=[up], colors='b', alpha=0.01)
    x, y = Globus_PET.bound(Shotn, time, rgr, zgr, ugr, up, share=False)
    return {'rb': x, 'zb': y}

def calcEqu(Shotn, time, Ipl, bet0, alf1, alf2):
    runRes = runGlobusPET(Shotn, time, Ipl, bet0, alf1, alf2)

    # если вдруг равновесие вообще не считается
    if runRes == -1:
        '''alf10 = alf1
        alf20 = alf2
        count = 0
        while runRes == -1 or count < 10:
            alf1, alf2 = round(random.uniform(alf10-0.2, alf10+0.2), 2), round(random.uniform(alf20-0.2, alf20+0.2), 2)
            runRes = runGlobusPET(Shotn, time, Ipl, bet0, alf1, alf2)
            count += 1
        if runRes == -1:'''
        return -1, -1, -1
    return alf1, alf2, runRes



def find_par3(Shotn, time, bet0, Ipl, alf1, alf2, plot=False, colorloc='g'):
    '''func for finding beta, alf1 and alf2'''

    Path_res = 'c:/work/Code/EqulibriumCodes/results/'

    MCC = MCCbound(Shotn, time)
    '''with open(Path_res + str(Shotn) + '_' + str(round((time / 1e3), 3)) + '_bound' + '.txt', 'r') as file2:
        rb = []
        zb = []
        for line in file2:
            data = line.split()
            rb.append(float(data[0]))
            zb.append((float(data[1])))


    PET = {'rb': rb, 'zb': zb}'''

    alf1, alf2, res = calcEqu(Shotn, time, Ipl, bet0, alf1, alf2)
    #print('result alpha:', alf1, alf2)

    if alf1 == -1:
        return -1

    PET = getPetBound(Shotn, time)
    chi2 = chi2_bound(MCC, PET)
    print(chi2, alf1, alf2, res)


    if plot:
        plt.plot(PET['rb'], PET['zb'], color=colorloc, label=r'$\alpha_1$=%.6f, $\alpha_2$=%.6f' %(alf1, alf2))
        plt.scatter(MCC['rb'], MCC['zb'], color='r')
    return chi2


Shotn, time = 44330, 190
bet0, Ipl = 0.2960, 301
psidia = 5.0932

alf1, alf2 = 1, 1
chi2 = find_par3(Shotn, time, bet0, Ipl, alf1, alf2, plot=False, colorloc='g')

petb = getPetBound(Shotn, time)
Mcc = MCCbound(Shotn, time)

plt.plot(petb['rb'], petb['zb'])
plt.plot(Mcc['rb'], Mcc['zb'])
plt.xlim(0, 1)
plt.ylim(-0.5, 0.5)
plt.show()

alf1, alf2 = 1, 1
for beta in [0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34]:
    chi2 = find_par3(Shotn, time, beta, Ipl, alf1, alf2, plot=False, colorloc='g')
    qdata, psi = Globus_PET.q_data_read(Shotn, time, PATH)
    chiPsi = (psidia - psi * 1e3) / psidia
    print('Result: ', beta, chi2, chiPsi)

stop

for alf2 in [0.7, 1, 1.2, 1.5]:
    for alf1 in [1.2, 1.4]:
        chi2 = find_par3(Shotn, time, bet0, Ipl, alf1, alf2, plot=False, colorloc='g')
        qdata, psi = Globus_PET.q_data_read(Shotn, time, PATH)
        chiPsi = (psidia - psi*1e3) / psidia
        print('Result: ', alf1, alf2, chi2, chiPsi)

'''Shotn, time = 44338, 170
allqata = []
alf2 = 1
plt.figure(figsize=(5, 8))
plt.xlim(0, 1)
plt.ylim(-0.8, 0.8)
alf1, alf2 = -0.03749999999999998, 0.9375
chi2 = find_par3(44338, 170, 0.2575, 297, alf1, alf2, plot=True, colorloc='g')
allqata.append(Globus_PET.q_data_read(Shotn, time, PATH))
alf1, alf2 = 0.125, 1.03125
find_par3(44338, 170, 0.2575, 297, alf1, alf2, plot=True, colorloc='b')
allqata.append(Globus_PET.q_data_read(Shotn, time, PATH))
alf1, alf2 = 0.56875, 1.375
find_par3(44338, 170, 0.2575, 297, alf1, alf2, plot=True, colorloc='m')
allqata.append(Globus_PET.q_data_read(Shotn, time, PATH))
alf1, alf2 = 0.66875, 1.5
find_par3(44338, 170, 0.2575, 297, alf1, alf2, plot=True, colorloc='c')
allqata.append(Globus_PET.q_data_read(Shotn, time, PATH))
alf1, alf2 = 0.8, 1.5
find_par3(44338, 170, 0.2575, 297, alf1, alf2, plot=True, colorloc='orange')
allqata.append(Globus_PET.q_data_read(Shotn, time, PATH))
plt.legend()

plt.figure()
for i, el in enumerate(allqata):
    plt.plot(el[0]['Psi'], el[0]['P'], label=i)
    print(el[1])
plt.legend()

plt.show()

#find_par3(44330, 168)

a1 = [0.4, 0.5, 0.6, 0.7, 0.8]
a2 = [0.5, 0.75, 1, 1.25, 1.5]
resa1 = []
resa2 = []
reschi = []

for i in range(5):
    alf1, alf2 = a1[i], a2[i]
    dalf0 = 0.5
    dalf = dalf0
    chi2res = 10
    eps = 0.044
    count = 0
    while chi2res > 0.03 and dalf >= 0.0625:
        print(chi2res, dalf)
        dalf = dalf/2

        alf2List = [alf2-dalf, alf2, alf2+dalf]
        chi2List = []
        for alf2 in alf2List:
            chi2List.append(abs(find_par3(44338, 170, 0.2575, 297, alf1, alf2)))
        chi2res = min(chi2List)
        alf2 = alf2List[chi2List.index(chi2res)]

        alf1List = [alf1 - dalf, alf1, alf1 + dalf]
        chi2List = []
        for alf1 in alf1List:
            chi2List.append(abs(find_par3(44338, 170, 0.2575, 297, alf1, alf2)))
        chi2res = min(chi2List)
        alf1 = alf1List[chi2List.index(chi2res)]

        count+=1
        if chi2res == 1:
            continue
        print('step:', count, chi2res, alf1, alf2)
    reschi.append(chi2res)
    resa1.append(alf1)
    resa2.append(alf2)
for i in range(5):
    print('RESULT: ', reschi[i], resa1[i], resa2[i])

plt.show()'''