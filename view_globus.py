import matplotlib.pyplot as plt
import json
from scipy import interpolate as interp
import numpy as np

PATH = 'c:/work/equilibrium/Globus_PET/'

def compare_bound():
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

    rgr = [float(i) for i in out_data[7:ni+7]]
    zgr = [float(i) for i in out_data[ni+7:ni+nj+7]]

    nugr = ni + nj + 7

    ugr = []
    curf = []

    for j in range(nj):
        ugr.append(out_data[nugr+ni*(j):nugr+ni*(j+1)])


    ncurf=nugr+2*ni*nj
    nipr=ncurf+ni*nj


    rm=out_data[nipr]
    zm=out_data[nipr+1]
    um=out_data[nipr+2]
    rx0=out_data[nipr+3]
    zx0=out_data[nipr+4]
    ux0=out_data[nipr+5]
    up=float(out_data[nipr+6])


    rxb=out_data[nipr+7:nipr+7+nxb]
    zxb=out_data[nipr+7+nxb:nipr+7+2*nxb]

    #slen=sum(((rxb(2:nxb)-rxb(1:nxb-1)).^2+(zxb(2:nxb)-zxb(1:nxb-1)).^2).^0.5)

    plt.figure(figsize=(5,8))
    cs = plt.contour(rgr, zgr, ugr, levels=[up], colors='b')
    plt.xlim(0,1)
    plt.ylim(-0.7, 0.7)
    plt.grid()
    bound = cs.allsegs
    x = []
    y = []
    for i in bound[0][0]:
        #print(i)
        x.append(float(i[0]))
        y.append(float(i[1]))
    #plt.figure(figsize=(5,8))
    plt.plot(x, y)
    '''plt.xlim(0,1)
    plt.ylim(-0.7, 0.7)
    plt.grid()'''


    with open(PATH + 'MCC.json', 'r') as file2:
        mcc_bound = json.load(file2)

    plt.plot([i/100 for i in mcc_bound['r']], [i/100 for i in mcc_bound['z']], 'r')

    diff_x = np.diff(x)
    diff_y = np.diff(y)
    area=abs(0.5*sum([y[i]*diff_x[i] - x[i]*diff_y[i] for i in range(len(x) - 1)]))

    diff_r = np.diff([i/100 for i in mcc_bound['r']])
    diff_z = np.diff([i/100 for i in mcc_bound['z']])
    area2 = abs(0.5*sum([[i/100 for i in mcc_bound['z']][j]*diff_r[j] - [i/100 for i in mcc_bound['r']][j]*diff_z[j] for j in range(len(mcc_bound['r']) - 1)]))

    print('-------------------')
    print(area, area2)
    print('-------------------')
    plt.show()


