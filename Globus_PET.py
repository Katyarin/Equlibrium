import ripper
import Find_boundary
import os
import niifaRipper
import numpy as np
import time as timer

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

#timer.sleep(10)
#os.system('"c:/work/equilibrium/Globus_PET/helisvf.exe"')







