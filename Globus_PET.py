import ripper
import Find_boundary
import os
import niifaRipper
import numpy as np
import time as timer
import matplotlib.pyplot as plt
import json
import subprocess

PATH = 'c:/work/equilibrium/Globus_PET/'

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def get_coils(Shotn, time):
    signals = niifaRipper.extract_niifa('//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/magn_data/', Shotn)

    i_need = [ 'Ipf1', 'Ipf3', 'Ics', 'Ipl', 'Icc', 'Ihfc', 'Ivfc']

    signals['Ipf2'] = [0] * len(signals['t_ms'])

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


def COIL_upd(Shotn, time, I_coil):
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

def compare_bound(par, ax, show=False, compare='dot', inside=False, share=False):
    out_data = []
    with open(PATH + 'out.wr', 'r') as file:
        for line in file:
            out_data.extend(line.split())

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

    for j in range(nj):
        ugr.append(out_data[nugr + ni * j:nugr + ni * (j + 1)])
        curf.append(out_data[nugr + ni * nj + ni * j: nugr + ni*nj + ni*(j + 1)])

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
        ax.set_ylim(-0.7, 0.7)
        ax.grid()
    cs = ax.contour(rgr, zgr, ugr, levels=[up], colors='b', alpha=0.1)
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
    for i in bound[0][bound_ind]:
        # print(i)
        x.append(float(i[0]))
        y.append(float(i[1]))
    # plt.figure(figsize=(5,8))
    if show == True:
        ax.plot(x, y, 'g')

    with open('bound.txt', 'w') as bndfile:
        for i, ri in enumerate(x):
            bndfile.write(' %8.5f ' % ri)
            bndfile.write(' %8.5f ' % y[i])
            bndfile.write('\n')
    '''plt.xlim(0,1)
    plt.ylim(-0.7, 0.7)
    plt.grid()'''

    with open(PATH + 'MCC.json', 'r') as file2:
        mcc_bound = json.load(file2)

    if show == True:
        ax.plot([i / 100 for i in mcc_bound['r']], [i / 100 for i in mcc_bound['z']], 'r')

    with open(PATH + 'BLANFW.DAT', 'r') as file3:
        blanfw = []
        for line in file3:
            blanfw.append(line.split())

    Rvv = float(blanfw[0][0])
    nvv = int(blanfw[1][0])

    ax.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
    ax.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')

    with open(PATH + 'LIMPNT.dat', 'r') as file4:
        limpnt = []
        for line in file4:
            limpnt.append(line.split())

    nlim = int(limpnt[0][0])

    ax.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
             [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')

    ax.scatter(rx0, zx0, marker='x')
    ax.scatter(rm, zm, marker='x')

    with open('dots.txt', 'w') as dotfile:
        dotfile.write('magnetic axis    ')
        dotfile.write(' %8.5f ' % rm)
        dotfile.write(' %8.5f ' % zm)
        dotfile.write('\n')
        dotfile.write('x-dot    ')
        dotfile.write(' %8.5f ' % rx0)
        dotfile.write(' %8.5f ' % zx0)


    print('_________compare bounds______________________________')
    print('center: ', rm, zm)
    print('x-dot: ', rx0, zx0)
    print('x_max = ', max(x))
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
        return area2 - area, (abs(dif_y1) + abs(dif_y2)) / 2
    elif compare == 'dot':
        dif_x1 = abs(max(x)) - abs(max([i / 100 for i in mcc_bound['r']]))
        dif_x2 = abs(min(x)) - abs(min([i / 100 for i in mcc_bound['r']]))
        dif_y1 = abs(max(y)) - abs(max([i / 100 for i in mcc_bound['z']]))
        dif_y2 = abs(min(y)) - abs(min([i / 100 for i in mcc_bound['z']]))
        print('-------------------')
        print((abs(dif_y1) + abs(dif_y2)) / 2)
        print(dif_x1, dif_x2)
        print('-------------------')
        return (abs(dif_x1) + abs(dif_x2)) / 2, (abs(dif_y1)) / 2
    else:
        print('ERROR')
        return 0, 0

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

def find_par(par, Shotn, time, I_coil, betta_po, li, bounds, show2=False):
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
            dif_x, dif_y = compare_bound(change / 100, show=show2)
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


def find_par2(par, Shotn, time, I_coil, betta_po, li, bounds, pdf, show2=False, how='not_all'):
    dif_list = []
    dif_list2 = []
    dif_list3 = []
    min_par = 0
    min_par_want = 0
    minimum = 1000

    res = {'li': [], 'bp': []}
    alf11, alf22 = find_li2(li)
    k = 0
    fig = plt.figure()
    fig.suptitle('t = ' + str(time) + ' bI = ' + str(betta_po))
    pic = 0
    fig.set_figheight(10)
    fig.set_figwidth(20)
    for change in range(bounds[0], bounds[1], bounds[2]):
        if par == 'betta_po':
            betta_po = change / 100
            print(betta_po)
        elif par == 'li':
            alf11, alf22 = find_li2(change / 100)
            print('li_want ', change / 100)
            print(alf11, alf22)
        else:
            print('error!', par)
            break
        DURS_upd(Shotn, time, I_coil['Ipl'], betta_po, alf11, alf22)

        try:
            process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate(timeout=5)
            result = EVOL_res()

            pic += 1
            ax = fig.add_subplot(2, 7, pic)

            if result == -1:
                print('NOT COUNT')
                continue
            if par == 'betta_po':
                print('res: ', result['betpol'])
                name = result['betpol']
            elif par == 'li':
                print('Li_res: ', result['li'])
                name = 'li = ' + str(round(result['li'], 2)) + ', bp = ' + str(round(result['betpol'], 2))
            #ax.set_title(str(name))
            dif_x, dif_y = compare_bound(name, ax, show=show2)
            if par == 'betta_po':
                dif = dif_y
            else:
                dif = (abs(dif_x)+abs(dif_y))/2
            dif_list.append(abs(dif_x))
            dif_list2.append((dif_y))
            dif_list3.append(dif)
            res['li'].append(result['li'])
            res['bp'].append(result['betpol'])
            if minimum > abs(dif):
                minimum = abs(dif)
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
    print('for %s min dif value %f with par value %f' % (par, minimum, min_par))
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
        axes[2].scatter(min_par,  minimum, color='r')
        pdf.savefig(fig2)

    return min_par, min_par_want, res


def find_bound(Shotn, time, I_coil, betta_I, li, alf1=0, k=0,  show2=True, inside=False):
    alf11, alf22 = find_li2(li, alf1=alf1)
    print("'''''''''alpha''''''''''''''")
    print(alf11, alf22)
    #print(Shotn, time, I_coil, betta_I, alf11, alf22)
    DURS_upd(Shotn, time, I_coil['Ipl'], betta_I, alf11, alf22)
    try:
        process = subprocess.Popen(["run.bat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(timeout=5)
        result = EVOL_res()
        if result == -1:
            print('NOT COUNT')
        else:
            name = '#' + str(Shotn) + ', time = ' + str(time) + ', bI = ' + str(betta_I) + '\n' + 'li = ' + str(round(result['li'], 2)) + ', bp = ' + str(round(result['betpol'], 2))
            fig = plt.figure(figsize=(6,10))
            ax = fig.add_subplot(1,1,1)
            dif_x, dif_y = compare_bound(name, ax, show=show2, inside=inside)
            plt.savefig('plots/' + str(Shotn) + '_' + str(li) + '_' + str(betta_I) + '_' + str(k) + '.png', dpi=600)

    except subprocess.TimeoutExpired:
        print('time over')
        subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))




