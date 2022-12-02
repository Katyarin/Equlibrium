import shtRipper
import Find_boundary
import os
import niifaRipper
import numpy as np
import time as timer
import matplotlib.pyplot as plt
import json
import subprocess
from scipy import interpolate as inter
import get_TS
from pathlib import Path
pi = 3.14159265359

PATH = 'c:/work/equilibrium/Globus_PET/'
sht_PATH = '//172.16.12.127/Data/sht'
PATH_res = 'results/'

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def get_coils(Shotn, time, shotn_sub=0):
    signals = niifaRipper.extract_niifa('//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/magn_data/', Shotn)
    print(signals.keys())



    i_need = ['Ipf1', 'Ipf3', 'Ics', 'Ipl', 'Icc', 'Ihfc', 'Ivfc']

    signals['Ipf2'] = [0] * len(signals['t_ms'])

    file_curr = Path('currents/' + str(Shotn) + '.txt')
    if not file_curr.is_file():
        with file_curr.open(mode='w') as file:
            file.write('%10s' %'time')
            for i in i_need:
                file.write('%10s' %i)
            file.write('\n')

    I_coil = {}
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
            I_coil[I] = smooth(signals[I], 95)[index_time2] #kA
        else:
            I_coil[I] = signals[I][index_time2] #kA
    with file_curr.open(mode='a') as file:
        file.write('%10.4f' % time)
        for i in i_need:
            file.write('%10.4f' %I_coil[i])
        file.write('\n')
    return I_coil


def DURS_upd(Shotn, time, Ip, betta_po, alf11, alf22):
    DURS = []
    with open(PATH + 'DURS.DAT', 'r') as file1:
        for line in file1:
            DURS.append(line.split())

    DURS_new = DURS
    DURS_new[0][3] = '#' + str(Shotn)
    DURS_new[0][4] = ' t=' + str(time) + 's,'
    DURS_new[0][5] = ' Ip=' + str(round(Ip / 1e3, 3))
    DURS_new[1][0] = str(betta_po) + 'd0'
    DURS_new[2][0] = str(round(Ip / 1e3, 3)) + 'd0'
    DURS_new[7][0] = str(alf11)
    DURS_new[8][0] = str(alf22)
    DURS_new[10][0] = str(alf11)
    DURS_new[11][0] = str(alf22)

    with open(PATH + 'DURS.DAT', 'w') as file2:
        for line in DURS_new:
            for element in range(len(line)):
                if element == 0:
                    file2.write(' %-12s' % line[element])
                else:
                    file2.write(' %s' % line[element])
            file2.write('\n')

    print('DURS.DAT updated')


def COIL_upd(Shotn, time, I_coil, Bt):
    COIL = []
    with open(PATH + 'COIL.DAT', 'r') as file3:
        for line in file3:
            COIL.append(line.split())

    COIL_new = COIL
    COIL_comments = []

    i = 0
    for line in COIL_new:
        COIL_comments.append([])
        for element in range(len(line)):
            if line[element] == '!' or line[element][0] == '!':
                COIL_comments[i]= line[element:]
                del line[element:]
                break
        i += 1

    for line in range(len(COIL_new)):
        if COIL_new[line] != []:
            COIL_new[line] = [float(s) for s in COIL_new[line]]

    COIL_comments[0][6] = '#' + str(Shotn)
    COIL_comments[0][7] = ' t=' + str(time) + 's,'
    COIL_comments[0][8] = ' Ip=' + str(round(I_coil['Ipl'] / 1e3, 3))


    for i in [3, 5]:
        COIL_new[i][6] = I_coil['Ipf1'] * 1e-3

    for i in [13, 15]:
        COIL_new[i][6] = I_coil['Ipf3'] * 1e-3

    for i in [18, 20, 22, 24]:
        COIL_new[i][6] = I_coil['Ihfc'] * 1e-3

    for i in [27, 29]:
        COIL_new[i][6] = I_coil['Ivfc'] * 1e-3

    for i in range(32,43,2):
        COIL_new[i][6] = I_coil['Icc'] * 1e-3

    for i in range(45,56,2):
        COIL_new[i][6] = I_coil['Ics'] * 1e-3

    with open(PATH + 'COIL.DAT', 'w') as file4:
        for line in range(len(COIL_new)):
            for element in range(len(COIL_new[line])):
                #print(el)
                el = COIL_new[line][element]
                if (abs(el) > 10 and len(COIL_new[line]) < 2) or element == 8:
                    file4.write(' %i' % el)
                elif abs(el) > 0.0099 and abs(el) < 1 and abs(el) != 0:
                    file4.write(' %8.5f' % el)
                elif abs(el) == 0 or abs(el) > 1:
                    file4.write(' %8.1f' % el)
                else:
                    file4.write(' %16.6e' % el)
            for el2 in COIL_comments[line]:
                #print(el2)
                file4.write(' %s' % el2)
            file4.write('\n')
        print('COIL.DAT updated')

    with open(PATH + 'DINA_ADD.DAT', 'w') as file5:
        file5.write(' RSO  ( cm)   BT0 (kGauss) ')
        file5.write('\n')
        file5.write(str(' 36.          ' + str(int(Bt*10)) + '.'))
        file5.write('\n')
        file5.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('DINA_ADD.DAT updated')


def EVOL_res():
    EVOL0 = []
    with open(PATH + 'EVOL0.PRT', 'r') as file2:
        for line in file2:
            EVOL0.append(line.split())
    if len(EVOL0) < 376:
        return -1
    else:
        Rc = (float(EVOL0[388][EVOL0[388].index('Rmax') + 2]) + float(EVOL0[389][EVOL0[389].index('Rmin') + 2])) / 2
        Result = {'li': float(EVOL0[377][EVOL0[377].index('LI3') + 2]), 'betpol': float(EVOL0[376][EVOL0[376].index('BETpol') + 2]),
                  'Rc': Rc}
        return Result

def compare_bound(Shotn, time, par, ax, show=False, curr=False, compare='dot', inside=False, share=False):

    out_data_raw = []

    with open(PATH + 'out.wr', 'r') as file:
        for line in file:
            out_data_raw.extend(line.split())


    q_data = {}
    line_count = 0
    with open(PATH + 'q.pr', 'r') as q_file:
        for line in q_file:
            data = line.split()
            if line_count < 1:
                line_count+=1
                continue
            elif line_count == 1:
                for i in data:
                    q_data[i] = []
                line_count += 1
            elif line_count < 63:
                i = 0
                for key in q_data.keys():
                    q_data[key].append(float(data[i]))
                    i+=1
                line_count += 1
            elif line_count == 65:
                Ftor_pl = float(data[4])
                line_count += 1
            else:
                line_count += 1

    #print(q_data.keys())
    P_axis = q_data['P'][0]
    out_data = out_data_raw
    for i in range(len(out_data_raw)):
        out_data[i] = float(out_data_raw[i])

    ni = int(out_data[0])
    nj = int(out_data[1])
    ni1 = int(out_data[2])
    nj1 = int(out_data[3])
    ni2 = int(out_data[4])
    nj2 = int(out_data[5])
    nxb = int(out_data[6])

    rgr = [float(i) for i in out_data[7:ni + 7]]
    zgr = [float(i) for i in out_data[ni + 7:ni + nj + 7]]

    nugr = ni + nj + 7

    ugr = []
    curf = []
    ipr = []

    for j in range(nj):
        ugr.append(out_data[nugr + ni * j:nugr + ni * (j + 1)])
        curf.append(out_data[nugr + ni * nj + ni * j: nugr + ni * nj + ni * (j + 1)])
        ipr.append(out_data[nugr + 2 * ni * nj + ni *j: nugr + 2 * ni * nj + ni * (j + 1)])

    norm_ugr = []
    for i, elem in enumerate(ugr):
        norm_ugr.append([])
        for j in elem:
            norm_ugr[i].append(j * 8 * pi * pi / 10)

    if share:
        with open(PATH_res + str(Shotn) + '_' + str(round(time, 3)) + '_eq_res.txt', 'w') as resfile:
            resfile.write('          ')
            for i, z in enumerate(zgr):
                resfile.write(' %8.5f ' % z)
            resfile.write('\n')
            for i, r in enumerate(rgr):
                resfile.write(' %8.5f ' % r)
                for j in range(len(zgr)):
                    resfile.write(' %8.5f ' % norm_ugr[i][j])
                resfile.write('\n')

        with open(PATH_res + str(Shotn) + '_' + str(round(time, 3)) + '_ipr_res.txt', 'w') as resfile:
            resfile.write('          ')
            for i, z in enumerate(zgr):
                resfile.write(' %8.5f ' % z)
            resfile.write('\n')
            for i, r in enumerate(rgr):
                resfile.write(' %8.5f ' % r)
                for j in range(len(zgr)):
                    resfile.write(' %8.5f ' % ipr[i][j])
                resfile.write('\n')

    ncurf = nugr + 2 * ni * nj
    nipr = ncurf + ni * nj

    rm = float(out_data[nipr])
    zm = float(out_data[nipr + 1])
    um = out_data[nipr + 2]
    rx0 = float(out_data[nipr + 3])
    zx0 = float(out_data[nipr + 4])
    ux0 = out_data[nipr + 5]
    up = float(out_data[nipr + 6])

    rxb = out_data[nipr + 7:nipr + 7 + nxb]
    zxb = out_data[nipr + 7 + nxb:nipr + 7 + 2 * nxb]

    # slen=sum(((rxb(2:nxb)-rxb(1:nxb-1)).^2+(zxb(2:nxb)-zxb(1:nxb-1)).^2).^0.5)

    if show == True:
        #plt.figure(figsize=(5, 8))
        ax.set_title(par)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.8, 0.8)
        ax.grid()
    if show:
        cs = ax.contour(rgr, zgr, ugr, levels=[up], colors='b', alpha=0.01)
    else:
        cs = plt.contour(rgr, zgr, ugr, levels=[up], colors='b', alpha=0.01)
    if inside:
        ax.contour(rgr, zgr, ugr, levels=20, colors='b', alpha=0.1)
    bound = cs.allsegs
    x = []
    y = []
    bound_ind = 0
    if len(bound[0]) > 1:
        max_bound = 0
        index_max = 0
        for i in range(len(bound[0])):
            if max_bound < len(bound[0][i]):
                max_bound = len(bound[0][i])
                index_max = i
        bound_ind = index_max
    print(bound_ind)
    for i in bound[0][bound_ind]:
        # print(i)
        x.append(float(i[0]))
        y.append(float(i[1]))
    # plt.figure(figsize=(5,8))
    if show == True:
        ax.plot(x, y, 'g', label='PET')




    if share:
        with open(PATH_res + str(Shotn) + '_' + str(round(time, 3)) + '_bound.txt', 'w') as bndfile:
            for i, ri in enumerate(x):
                bndfile.write(' %8.5f ' % ri)
                bndfile.write(' %8.5f ' % y[i])
                bndfile.write('\n')

    with open(PATH + 'MCC.json', 'r') as file2:
        mcc_bound = json.load(file2)

    if show == True:
        ax.plot([i / 100 for i in mcc_bound['r']], [i / 100 for i in mcc_bound['z']], 'r', label='MCC')
        ax.legend()
    with open(PATH + 'BLANFW.DAT', 'r') as file3:
        blanfw = []
        for line in file3:
            blanfw.append(line.split())

    Rvv = float(blanfw[0][0])
    nvv = int(blanfw[1][0])

    if show:

        ax.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
        ax.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')

    with open(PATH + 'LIMPNT.dat', 'r') as file4:
        limpnt = []
        for line in file4:
            limpnt.append(line.split())

    nlim = int(limpnt[0][0])
    r_lim = [float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])]
    z_lim = [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])]

    strike_point = {'inner': [], 'outer': []}
    if min(y) < zx0:
        for j, ind in enumerate([len(x) - 1, 0]):
            for i in range(len(r_lim)):
                if r_lim[i] > x[ind] and r_lim[i + 1] < x[ind]:
                    k2, b2 = np.polyfit(r_lim[i:i + 2], z_lim[i:i + 2], 1)
            print(k2, b2)

            if j:
                ind_first = 0
                ind_last = ind + 6
            else:
                ind_first = ind - 5
                ind_last = ind + 1
            k1, b1 = np.polyfit(x[ind_first:ind_last], y[ind_first:ind_last], 1)

            strike_point[list(strike_point.keys())[j]].extend([(b1 - b2) / (k2 - k1), (k2 * b1 - k1 * b2) / (k2 - k1)])
    else:
        for j in range(2):
            strike_point[list(strike_point.keys())[j]].extend([-1000,-1000])
    if show:
        ax.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
                 [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')

        ax.scatter(rx0, zx0, marker='x')
        ax.scatter(rm, zm, marker='x')
        for j in range(2):
            ax.scatter(strike_point[list(strike_point.keys())[j]][0], strike_point[list(strike_point.keys())[j]][1], marker='*')

    if share:
        with open(PATH_res + str(Shotn) + '_' + str(round(time, 3)) + '_dots.txt', 'w') as dotfile:
            dotfile.write('magnetic axis    ')
            dotfile.write(' %8.5f ' % rm)
            dotfile.write(' %8.5f ' % zm)
            dotfile.write('\n')
            dotfile.write('x-dot    ')
            dotfile.write(' %8.5f ' % rx0)
            dotfile.write(' %8.5f ' % zx0)

    """CURRENT"""
    zst = 1000
    for i in range(len(zgr)):
        if zst > abs(zm - zgr[i]):
            zst = abs(zm - zgr[i])
            jst = i

    #print(zst, jst)

    zst = 1000
    for i in range(len(zgr)):
        if zst > abs(0 - zgr[i]):
            zst = abs(0 - zgr[i])
            js0 = i

    #print(zst, js0)

    k = -1
    curf2 = []
    rsn = []
    for i in range(ni):
        if abs(float(curf[jst][i])) > 1e-5:
            k = k + 1
            if k == 0:
                curf2.append(float(curf[jst][i - 1]))
                rsn.append(rgr[i - 1])
            curf2.append(float(curf[jst][i]))
            rsn.append(rgr[i])
            ki = i

    curf2.append(float(curf[jst][ki + 1]))
    rsn.append(rgr[ki + 1])

    sum1 = []
    for i in range(len(curf)):
        sum1.append(sum(curf[i]))

    Ipl = sum(sum1) * (rgr[1] - rgr[0]) * (zgr[1] - zgr[0])

    if curr:
        plt.figure()
        plt.plot(rsn, [i / Ipl for i in curf2])
        plt.grid()
        plt.ylim(0, 15)
        plt.xlim(0.1, 0.7)
        plt.ylabel('MA/(m*m)')
        plt.xlabel('r(m)')

    """PRESSURE"""
    P_psi = inter.interp1d(q_data['Psi'], q_data['P'], kind='cubic')

    press = []

    for i, psi in enumerate(norm_ugr):
        press.append([])
        for j in psi:
            if max(q_data['Psi']) > j > min(q_data['Psi']):
                press[i].append(P_psi(j))
            elif max(q_data['Psi']) < j:
                press[i].append(P_psi(max(q_data['Psi'])))
            else:
                press[i].append(0)

    W_prob = []
    V_prob = []

    for z in range(len(press[0])):
        W_prob.append(np.trapz([press[r][z] * 1e6 for r in range(len(press))], zgr))
    W_all = 3 / 2 * np.trapz([W * 2 * pi * rgr[r] for r, W in enumerate(W_prob)], rgr)

    for z in range(len(ipr[0])):
        V_prob.append(np.trapz([ipr[r][z] for r in range(len(ipr))], zgr))
    V = np.trapz([v * 2 * pi * rgr[r] for r, v in enumerate(V_prob)], rgr)
    S = np.trapz([v for r, v in enumerate(V_prob)], rgr)

    """q"""
    q_psi = inter.interp1d(q_data['Psi'], q_data['qplr'], kind='cubic')

    q_q = []

    q_95 = q_psi(up* 8 * pi * pi / 10 *0.95)


    for i, psi in enumerate(norm_ugr):
        q_q.append([])
        for j in psi:
            if max(q_data['Psi']) > j > min(q_data['Psi']):
                q_q[i].append(q_psi(j))
            elif max(q_data['Psi']) < j:
                q_q[i].append(q_psi(max(q_data['Psi'])))
            else:
                q_q[i].append(0)

    if share:
        with open(PATH_res + str(Shotn) + '_' + str(round(time, 3)) + '_q_res.txt', 'w') as qresfile:
            qresfile.write('          ')
            for i, z in enumerate(zgr):
                qresfile.write(' %8.5f ' % z)
            qresfile.write('\n')
            for i, r in enumerate(rgr):
                qresfile.write(' %8.5f ' % r)
                for j in range(len(zgr)):
                    qresfile.write(' %8.5f ' % q_q[i][j])
                qresfile.write('\n')




    '''print('_________compare bounds______________________________')
    print('center: ', rm, zm)
    print('x-dot: ', rx0, zx0)
    print('x_max = ', max(x))'''
    if compare == 'area':
        diff_x = np.diff(x)
        diff_y = np.diff(y)
        area = abs(0.5 * sum([y[i] * diff_x[i] - x[i] * diff_y[i] for i in range(len(x) - 1)]))

        diff_r = np.diff([i / 100 for i in mcc_bound['r']])
        diff_z = np.diff([i / 100 for i in mcc_bound['z']])
        area2 = abs(0.5 * sum(
            [[i / 100 for i in mcc_bound['z']][j] * diff_r[j] - [i / 100 for i in mcc_bound['r']][j] * diff_z[j] for j in
            range(len(mcc_bound['r']) - 1)]))

        print('-------------------')
        print(area, area2)
        print('-------------------')
        dif_y1 = abs(max(y)) - abs(max([i / 100 for i in mcc_bound['z']]))
        dif_y2 = abs(min(y)) - abs(min([i / 100 for i in mcc_bound['z']]))
        '''if show == True:
            plt.show()'''
        return area2 - area, (abs(dif_y1) + abs(dif_y2)) / 2, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q_95
    elif compare == 'dot':
        dif_x1 = abs(max(x)) - abs(max([i / 100 for i in mcc_bound['r']]))
        dif_x2 = abs(min(x)) - abs(min([i / 100 for i in mcc_bound['r']]))
        dif_y1 = abs(max(y)) - abs(max([i / 100 for i in mcc_bound['z']]))
        dif_y2 = abs(min(y)) - abs(min([i / 100 for i in mcc_bound['z']]))
        '''print('-------------------')
        print((abs(dif_y1) + abs(dif_y2)) / 2)
        print(dif_x1, dif_x2)
        print('-------------------')'''
        return (abs(dif_x1) + abs(dif_x2)) / 2, (abs(dif_y1)) / 2, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q_95
    else:
        print('ERROR')
        return 0, 0, [], [], [], 0, 0, 0, 0, 0, {}, 0

def find_li(li):
    li_arr = np.loadtxt('li_res.txt')
    al11_arr = [i / 100 for i in range(60, 120)]
    al22_arr = [i / 100 for i in range(80, 150)]
    ind_find = np.where(np.logical_and(li_arr >= li - 0.002, li_arr <= li + 0.002))
    if max(ind_find[0]) < 30:
        max_i = max(ind_find[0])
        ind_0 = np.where(ind_find[0] == max_i)[0]
    else:
        ind = np.where((ind_find[0] > 30) & (ind_find[0] < 50))
        ind_0 = ind[0][int(len(ind[0]) / 2)]
        print(ind[0][int(len(ind[0]) / 2)])
    ind1 = ind_find[0][ind_0]
    ind2 = ind_find[1][ind_0]
    print('li: ', float(li_arr[ind1, ind2]))
    return al11_arr[int(ind1)], al22_arr[int(ind2)]

def find_li2(li, alf1=0):
    if alf1 == 0:
        alf1 = (li - 0.2 + 0.6) / 2
    else:
        alf1 = (li - 0.2 + 0.6)
    a1 = -0.049
    a2 = 1.34
    b11 = 0.41
    b12 = - 0.944
    b21 = -0.1345
    b22 = 0.26
    alf2 = (li - a1 - b11 * alf1 - b21 * alf1 * alf1) / (a2 + b12 * alf1 + b22 * alf1 * alf1)
    return alf1, alf2

def find_par(par, Shotn, time, I_coil, betta_po, li, bounds, show2=False, share=False):
    dif_list = []
    dif_list2 = []
    min_par = 0
    minimum = 1000
    res = {'li': [], 'bp': []}
    #alf11, alf22 = find_li(li)
    for change in range(bounds[0], bounds[1], bounds[2]):
        if par == 'betta_po':
            betta_po = change / 100
            print(betta_po)
        elif par == 'li':
            alf11, alf22 = find_li(change / 100)
            print('li_want ', change / 100)
            print(alf11, alf22)
        else:
            print('error!', par)
            break
        DURS_upd(Shotn, time, I_coil['Ipl'], betta_po, alf11, alf22)

        try:
            process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate(timeout=5)
            # print('out: ', stdout.decode('utf-8'))
            result = EVOL_res()
            if result == -1:
                print('NOT COUNT')
                continue
            dif_x, dif_y, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q95 = compare_bound(Shotn, time, change / 100, show=show2, share=share)
            if par == 'betta_po':
                dif = dif_y
            else:
                dif = dif_x
            dif_list.append(abs(dif_x))
            dif_list2.append((dif_y))
            if par == 'betta_po':
                print('res: ', result['betpol'])
            elif par == 'li':
                print('Li_res: ', result['li'])
            if minimum > abs(dif):
                minimum = abs(dif)
                if par == 'betta_po':
                    min_par = result['betpol']
                elif par == 'li':
                    min_par = result['li']
            res['li'].append(result['li'])
            res['bp'].append(result['betpol'])
        except subprocess.TimeoutExpired:
            print('time over')
            subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))
        '''print('err: ', stderr.decode('utf-8'))'''
    print('for %s min dif value %f with par value %f' %(par, minimum, min_par))


    if show2 == True:
        if par == 'betta_po':
            x = res['bp']
        elif par == 'li':
            x = res['li']
        else:
            x = [i / 100 for i in range(bounds[0], bounds[0] + len(dif_list))]

        plt.figure()
        plt.title('x difference')
        plt.plot(x, dif_list, 'o')
        plt.grid()

        plt.figure()
        plt.title('y difference')
        plt.plot(x, dif_list2, 'o')
        plt.grid()
        plt.show()

    return min_par, res


def find_par2(par, mode, Shotn, time, I_coil, betta_po, li, bounds, pdf, show2=False, how='not_all', psi_dia_sakh=0, psi_mode = 0):
    dif_list = []
    dif_list2 = []
    dif_list3 = []
    dif_list4 = []
    min_par = 0
    min_par_want = 0
    minimum = 1000

    res = {'li': [], 'bp': []}
    alf11, alf22 = find_li2(li)
    k = 0
    pic = 0
    if show2:
        fig = plt.figure()
        fig.suptitle('t = ' + str(time) + ' bI = ' + str(betta_po))
        fig.set_figheight(10)
        fig.set_figwidth(20)
    if mode == 'sep':
        minimum_delta = 1
        for change in range(bounds[0], bounds[1], bounds[2]):
            print(change / 100)
            if par == 'betta_po':
                betta_po = change / 100
                #print(betta_po)
            elif par == 'li':
                alf11, alf22 = find_li2(change / 100)
                '''print('li_want ', change / 100)
                print(alf11, alf22)'''
            else:
                print('error!', par)
                break
            DURS_upd(Shotn, time, I_coil['Ipl'], betta_po, alf11, alf22)

            try:
                process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate(timeout=5)
                result = EVOL_res()

                if show2:
                    pic += 1
                    ax = fig.add_subplot(2, 7, pic)
                else:
                    ax = 0

                if result == -1:
                    print('NOT COUNT')
                    continue
                if par == 'betta_po':
                    #print('res: ', result['betpol'])
                    name = result['betpol']
                elif par == 'li':
                    #print('Li_res: ', result['li'])
                    name = 'li = ' + str(round(result['li'], 2)) + ', bp = ' + str(round(result['betpol'], 2))
                #ax.set_title(str(name))
                dif_x, dif_y, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q95 = compare_bound(Shotn, time, name, ax, show=show2)
                if par == 'betta_po':
                    dif = dif_y
                else:
                    dif = (abs(dif_x)+abs(dif_y))/2
                dif_list.append(abs(dif_x))
                dif_list2.append((dif_y))
                dif_list3.append(dif)
                res['li'].append(result['li'])
                res['bp'].append(result['betpol'])
                psi_pet = Ftor_pl * 1000
                delta_psi = abs(psi_dia_sakh - psi_pet)
                if minimum > abs(dif):
                    minimum = abs(dif)
                    minimum_delta = delta_psi
                    if par == 'betta_po':
                        min_par = result['betpol']
                    elif par == 'li':
                        min_par = result['li']
                        min_par_want = change / 100
                else:
                    if how == 'all' or k < 2:
                        k += 1
                        continue
                    else:
                        break

            except subprocess.TimeoutExpired:
                print('time over')
                subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))
            '''print('err: ', stderr.decode('utf-8'))'''
    elif mode == 'psi':
        bnd_s = bounds[0]
        bnd_e = bounds[1]
        delta_psi = 1
        minimum_delta = 100
        ground_psi = 0.005
        if psi_mode:
            ground_psi = 0.05
        while delta_psi > ground_psi:
            par_test = round((bnd_e + bnd_s) / 2, 2)
            print(par_test)
            if par == 'betta_po':
                betta_po = par_test / 100
                #print(betta_po)
            elif par == 'li':
                alf11, alf22 = find_li2(par_test / 100)
                '''print('li_want ', par_test / 100)
                print(alf11, alf22)'''
            else:
                print('error!', par)
                break
            DURS_upd(Shotn, time, I_coil['Ipl'], betta_po, alf11, alf22)

            try:
                pic += 1
                print('.................................', betta_po, pic, '.................................')
                if pic < 11:
                    if show2:
                        ax = fig.add_subplot(3, 7, pic)
                    else:
                        ax = 0
                else:
                    print(delta_psi)
                    print('too many steps')
                    break



                process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate(timeout=5)
                result = EVOL_res()



                if result == -1:
                    #print('NOT COUNT')
                    continue
                if par == 'betta_po':
                    #print('res: ', result['betpol'])
                    name = result['betpol']
                elif par == 'li':
                    #print('Li_res: ', result['li'])
                    name = 'li = ' + str(round(result['li'], 2)) + ', bp = ' + str(round(result['betpol'], 2))
                #ax.set_title(str(name))
                dif_x, dif_y, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q95 = compare_bound(Shotn, time, name, ax, show=show2)
                if par == 'betta_po':
                    dif = dif_y
                else:
                    dif = (abs(dif_x)+abs(dif_y))/2
                dif_list.append(abs(dif_x))
                dif_list2.append((dif_y))
                dif_list3.append(dif)
                res['li'].append(result['li'])
                res['bp'].append(result['betpol'])

                psi_pet = Ftor_pl * 1000
                delta_psi = abs(psi_dia_sakh - psi_pet)
                dif_list4.append(delta_psi)
                '''if dif > 0.02:
                    delta_psi = 0.2'''
                '''print('!!!!!!!!!!!!!!!!!')
                print(psi_dia_sakh, psi_pet)
                print(delta_psi)
                print(bnd_s, bnd_e, par_test)
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')'''
                if delta_psi < minimum_delta:
                    minimum_delta = delta_psi
                    minimum = abs(dif)
                    if par == 'betta_po':
                        min_par = result['betpol']
                    elif par == 'li':
                        min_par = result['li']
                        min_par_want = par_test / 100
                if psi_dia_sakh > psi_pet:
                    bnd_s = par_test
                else:
                    bnd_e = par_test


            except subprocess.TimeoutExpired:
                print('time over')
                if par_test > 1.3 * 100:
                    bnd_e = par_test
                else:
                    bnd_s = par_test
                subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))
            '''print('err: ', stderr.decode('utf-8'))'''





    else:
        print('unknown mode!')
        stop
    #print('for %s min dif value %f with par value %f' % (par, minimum, min_par))
    if show2:
        pdf.savefig(fig)

    if show2 == True:
        if par == 'betta_po':
            x = res['bp']
        elif par == 'li':
            x = res['li']
        else:
            x = [i / 100 for i in range(bounds[0], bounds[0] + len(dif_list))]

        fig2, axes = plt.subplots(3, 1)
        axes[0].set_title('x difference')
        axes[0].plot(x, dif_list, 'o')
        axes[0].grid()

        axes[1].set_title('y difference')
        axes[1].plot(x, dif_list2, 'o')
        axes[1].grid()

        axes[2].set_title('av difference')
        axes[2].plot(x, dif_list3, 'o')
        axes[2].grid()
        axes[2].scatter(min_par,  minimum, color='r', marker='^', zorder=2.5)
        pdf.savefig(fig2)

        if mode == 'psi':
            fig3, axes2 = plt.subplots(1, 1)
            axes2.set_title('psi_dia - psi_pet, mV')
            axes2.plot(x, dif_list4, 'o')
            axes2.grid()
            pdf.savefig(fig3)

    return min_par, min_par_want, res, minimum_delta, minimum


def We(Shotn, time, rgr, zgr, norm_ugr, ipr, r_axis, Wi=False):
    q = 1.602176634e-19
    time = time * 1e3
    try:
        R, Te, ne_raw, P, error = get_TS.get_data(Shotn, time)  # tomson data
    except:
        return 0,0,0,0
    P_TS_raw = [i * 1e-6 for i in P['P']]
    R_TS_raw = [i / 1000 for i in R]
    ne = {'ne': [], 'err': []}

    ind_max_for_prof = len(R_TS_raw)

    for i, rr in enumerate(R_TS_raw):
        if rr > r_axis - 0.02:
            ind_max_for_prof = i+1

    P_TS = P_TS_raw[:ind_max_for_prof]
    R_TS = R_TS_raw[:ind_max_for_prof]
    ne['ne'] = ne_raw['ne'][:ind_max_for_prof]

    zst = 1000
    for i in range(len(zgr)):
        if zst > abs(0 - zgr[i]):
            zst = abs(0 - zgr[i])
            js0 = i
    Psi_R_for_TS = inter.interp1d(rgr, norm_ugr[js0], kind='cubic')
    Psi_R_for_TS_raw = inter.interp1d(rgr, norm_ugr[js0], kind='cubic')

    P_TS_psi = inter.interp1d(Psi_R_for_TS(R_TS), P_TS, kind='linear')

    ne_TS_psi = inter.interp1d(Psi_R_for_TS(R_TS), ne['ne'], kind='linear')

    ne_psi = []
    for i in norm_ugr[js0]:
        if max(Psi_R_for_TS(R_TS)) > i > min(Psi_R_for_TS(R_TS)):
            ne_psi.append(ne_TS_psi(i))
        elif i > max(Psi_R_for_TS(R_TS)):
            ne_psi.append(ne_TS_psi(max(Psi_R_for_TS(R_TS))))
        else:
            ne_psi.append(0)

    ne_r = inter.interp1d(rgr, ne_psi, kind='linear')

    ne_all_psi = []

    for i, psi in enumerate(norm_ugr):
        ne_all_psi.append([])
        for j in psi:
            if max(Psi_R_for_TS(R_TS)) > j > min(Psi_R_for_TS(R_TS)):
                ne_all_psi[i].append(ne_TS_psi(j))
            elif j > max(Psi_R_for_TS(R_TS)):
                ne_all_psi[i].append(ne_TS_psi(max(Psi_R_for_TS(R_TS))))
            else:
                ne_all_psi[i].append(0)

    TS_press = []


    for i, psi in enumerate(norm_ugr):
        TS_press.append([])
        for j in psi:
            if max(Psi_R_for_TS(R_TS)) > j > min(Psi_R_for_TS(R_TS)):
                TS_press[i].append(P_TS_psi(j))
            elif j > max(Psi_R_for_TS(R_TS)):
                TS_press[i].append(P_TS_psi(max(Psi_R_for_TS(R_TS))))
            else:
                TS_press[i].append(0)
    plt.figure()
    plt.plot(rgr, TS_press[js0], label='electron')
    plt.scatter(R_TS_raw, P_TS_raw, label='TS')
    plt.vlines(r_axis, min(P_TS_raw), max(P_TS_raw)*1.2, color='black')

    '''with open(str(Shotn) + '_' + str(time) + '_TS_press.json', 'w') as ts_file:
        for_wr = {'Pressure': [float(i) for i in TS_press[js0]], 'R': rgr, 'TS_data':
            {'R': R_TS, 'P_TS': P_TS, 'P_err': [i * 1e-6 for i in P['err']]}}
        json.dump(for_wr, ts_file)'''
    W_prob_TS = []
    ne_prob_TS = []

    V_prob = []
    for z in range(len(ipr[0])):
        V_prob.append(np.trapz([ipr[r][z] for r in range(len(ipr))], zgr))
    V = np.trapz([v * 2 * pi * rgr[r] for r, v in enumerate(V_prob)], rgr)

    for z in range(len(TS_press[0])):
        W_prob_TS.append(np.trapz([TS_press[r][z] * ipr[r][z]  * 1e6 for r in range(len(TS_press))], zgr))
        ne_prob_TS.append(np.trapz([ne_all_psi[r][z] * ipr[r][z] for r in range(len(TS_press))], zgr))
    W_electron = 3 / 2 * np.trapz([W * 2 * pi * rgr[r] for r, W in enumerate(W_prob_TS)], rgr)
    ne_av = np.trapz([n * 2 * pi * rgr[r] for r, n in enumerate(ne_prob_TS)], rgr) / V

    #print('average_ne: ', ne_av)
    if Wi:
        Z1 = 1
        Z2 = 6
        time_li = []
        Ri_li = {}
        Ti_li = {}
        Ti_err_li = {}

        try:
            l = 0
            with open('c:/work/equilibrium/Ti_data/' + str(Shotn) + '.txt', 'r') as Ti_file:
                for line in Ti_file:
                    data = line.split()
                    if data:
                        if l == 0:
                            time_li = [float(data[i]) for i in range(1, len(data), 2)]
                        elif l == 1:
                            t = 0
                            R_local = float(data[0])
                            for i in range(1, len(data), 2):
                                if data[i] != 'NaN':
                                    Ri_li[int(time_li[t] * 10)] = [R_local]
                                    Ti_li[int(time_li[t] * 10)] = [float(data[i])]
                                    Ti_err_li[int(time_li[t] * 10)] = [float(data[i+1])]
                                else:
                                    Ri_li[int(time_li[t] * 10)] = []
                                    Ti_li[int(time_li[t] * 10)] = []
                                    Ti_err_li[int(time_li[t] * 10)] = []
                                t+=1
                        else:
                            t = 0
                            #print(data)
                            R_local = float(data[0])
                            for i in range(1, len(data), 2):
                                if data[i] != 'NaN':
                                    Ri_li[int(time_li[t] * 10)].append(R_local)
                                    Ti_li[int(time_li[t] * 10)].append(float(data[i]))
                                    Ti_err_li[int(time_li[t] * 10)].append(float(data[i+1]))
                                t += 1
                        l+=1
                '''print(time_li)
                print(Ri_li)
                print(Ti_li)'''

            with open('c:/work/equilibrium/Zeff_data/' + str(Shotn) + '_Zeff.txt', 'r') as Zeff_file:
                Zeff_data = {'time': [], 'data': []}
                for line in Zeff_file:
                    data = line.split()
                    Zeff_data['time'].append(float(data[1]) * 1e3)
                    Zeff_data['data'].append(float(data[3]))
            #print(Zeff_data)
            if round(time, 1) in time_li:
                Ri_raw = Ri_li[int(round(time * 10))]
                Ti_raw = Ti_li[int(round(time * 10))]
                ind_max_for_prof_i = 0
                for i, rr in enumerate(Ri_raw):
                    if rr/1000 < r_axis - 0.010:
                        ind_max_for_prof_i = i+1
                Ri = Ri_raw[ind_max_for_prof_i:]
                Ti = Ti_raw[ind_max_for_prof_i:]
                # Ti_err = Ti_err_li[int(round(time * 10))]
                '''print(Ri, Ti)
                print('we are here')'''
                ni = ne_r([elem / 1000 for elem in Ri])
                ni_raw = ne_r([elem / 1000 for elem in Ri_raw])
                for i, el in enumerate(Zeff_data['time']):
                    if time - 2.5 < el < time + 2.5:
                        '''print(el)
                        print(Zeff_data['data'][i])'''
                        Zeff = Zeff_data['data'][i]
                #print('Zeff:', Zeff)
                ni_coeff = (Zeff - Z2) / (Z1 * Z1 - Z1 * Z2)
                #print(ni_coeff)
                Pi = [Ti[i] * ni[i] * ni_coeff * q * 1e-6 for i in range(len(Ri))]
                Pi_raw = [Ti_raw[i] * ni_raw[i] * ni_coeff * q * 1e-6 for i in range(len(Ri_raw))]
                # Pi_err = [Pi[i] * (Ti_err[i] / Ti[i]) for i in range(len(Ri))]

                Pi_psi = inter.interp1d(Psi_R_for_TS([elem / 1000 for elem in Ri]), Pi, kind='linear')
                ion_press = []
                for i, psi in enumerate(norm_ugr):
                    ion_press.append([])
                    for j in psi:
                        if max(Psi_R_for_TS([elem / 1000 for elem in Ri])) > j > min(
                                Psi_R_for_TS([elem / 1000 for elem in Ri])):
                            ion_press[i].append(Pi_psi(j))
                        elif j > max(Psi_R_for_TS([elem / 1000 for elem in Ri])):
                            ion_press[i].append(Pi_psi(max(Psi_R_for_TS([elem / 1000 for elem in Ri]))))
                        else:
                            ion_press[i].append(0)

                plt.plot(rgr, ion_press[js0], label='ion')
                plt.scatter([elem / 1000 for elem in Ri_raw], Pi_raw, label='CXRS')
                plt.legend()
                plt.grid()
                plt.ylim(0, max(Pi_raw)*1.2)
                plt.savefig('el_ion_prof/' + str(Shotn) + '_' + str(time) + '.png')

                W_prob_i = []
                for z in range(len(ion_press[0])):
                    W_prob_i.append(np.trapz([ion_press[r][z] * ipr[r][z]  * 1e6 for r in range(len(ion_press))], zgr))
                W_ion = 3 / 2 * np.trapz([W * 2 * pi * rgr[r] for r, W in enumerate(W_prob_i)], rgr)
            else:
                plt.legend()
                plt.grid()
                plt.ylim(0, max(P_TS_raw) * 1.2)
                plt.savefig('el_ion_prof/' + str(Shotn) + '_' + str(time) + '.png')
                W_ion = 0
                print('Ti file not found')
                W_ion = 0
        except FileNotFoundError:
            plt.legend()
            plt.grid()
            plt.ylim(0, max(P_TS_raw) * 1.2)
            plt.savefig('el_ion_prof/' + str(Shotn) + '_' + str(time) + '.png')
            W_ion = 0
            print('Ti file not found')

    else:
        plt.legend()
        plt.grid()
        plt.ylim(0, max(P_TS_raw) * 1.2)
        plt.savefig('el_ion_prof/' + str(Shotn) + '_' + str(time) + '.png')
        W_ion = 0

    return W_electron, error, W_ion, ne_av

def find_bound(Shotn, time, I_coil, betta_I, li, alf1=0, k=0, pdf=None, show2=True, inside=False, Wi=False, share=False):
    alf11, alf22 = find_li2(li, alf1=alf1)
    print("'''''''''alpha''''''''''''''")
    print(alf11, alf22)
    #print(Shotn, time, I_coil, betta_I, alf11, alf22)
    DURS_upd(Shotn, time, I_coil['Ipl'], betta_I, alf11, alf22)
    try:
        process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(timeout=5)
        result = EVOL_res()
        print(result)
        if result == -1:
            print('NOT COUNT')
        else:
            name = '#' + str(Shotn) + ', time = ' + str(time) + ', bI = ' + str(betta_I) + '\n' + 'li = ' + str(round(result['li'], 2)) + ', bp = ' + str(round(result['betpol'], 2))
            fig = plt.figure(figsize=(6,10))
            ax = fig.add_subplot(1,1,1)
            dif_x, dif_y, rgr, zgr, norm_ugr, W_all, V, S, P_axis, Ftor_pl, ipr, rm, strike_point, q95 = compare_bound(Shotn, time, name, ax, show=show2, inside=inside, share=share)
            We_el, error, W_ion, ne_av = We(Shotn, time, rgr, zgr, norm_ugr, ipr, rm, Wi)
            dif = (abs(dif_x) + abs(dif_y)) / 2
            if pdf:
                pdf.savefig(fig)
        return betta_I, str(round(result['li'], 3)), str(round(result['betpol'], 3)), W_all, We_el, V, S, P_axis, li, Ftor_pl, W_ion, dif, ne_av, strike_point, q95
    except subprocess.TimeoutExpired:
        print('time over')
        subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))




