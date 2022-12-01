import matplotlib.pyplot as plt
import numpy as np
result = []
eq_result = []
rb = []
zb = []
#PATH = 'c:/work/equilibrium/Globus_PET/'
Path_res = 'results/'
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 42567
time = 210 #ms

z_need = 0

plt.figure(figsize=(5, 8))
with open(PATH + 'BLANFW.DAT', 'r') as file3:
    blanfw = []
    for line in file3:
        blanfw.append(line.split())

Rvv = float(blanfw[0][0])
nvv = int(blanfw[1][0])

with open(PATH + 'LIMPNT.dat', 'r') as file4:
    limpnt = []
    for line in file4:
        limpnt.append(line.split())

with open(PATH + 'COIL.dat', 'r') as file7:
    coil_data = []
    for line in file7:
        coil_data.append(line.split())
r_coil = []
z_coil = []
dr_coil = []
dz_coil = []
coil_list = []
for line in coil_data:
    if len(line) == 11 or len(line) == 9 and line[0] != '!':
        r_coil.append(float(line[0]))
        z_coil.append(float(line[1]))
        dr_coil.append(float(line[2]))
        dz_coil.append(float(line[3]))
        #coil_list.append(line[10])
#print(coil_list)
nlim = int(limpnt[0][0])

"""with open(PATH + 'BLANFW.DAT', 'r') as file3:
    blanfw = []
    for line in file3:
        blanfw.append(line.split())

Rvv = float(blanfw[0][0])
nvv = int(blanfw[1][0])

plt.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
plt.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')

with open(PATH + 'LIMPNT.dat', 'r') as file4:
    limpnt = []
    for line in file4:
        limpnt.append(line.split())

nlim = int(limpnt[0][0])

plt.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
        [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')"""

with open(Path_res + str(Shotn) + '_' + str(time/1e3) + '_q_res.txt', 'r') as file:
    for line in file:
        result.append([float(i) for i in line.split()])

with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_eq_res.txt', 'r') as file:
    for line in file:
        eq_result.append([float(i) for i in line.split()])

with open(Path_res + str(Shotn) + '_' + str(time/1e3) + '_bound.txt', 'r') as file2:
    for line in file2:
        data = line.split()
        rb.append(float(data[0]))
        zb.append((float(data[1])))

dots = {}

with open(Path_res + str(Shotn) + '_' + str(time/1e3) + '_dots.txt', 'r') as file3:
    for line in file3:
        data = line.split()
        if len(data) < 4:
            dots[data[0]] = [float(data[1]), float(data[2])]
        else:
            dots[data[0] + '_' + data[1]] = [float(data[2]), float(data[3])]

zgr = result[0]
rgr = [result[i][0] for i in range(1, len(result))]
ugr = []
for i in range(len(result)):
    if i != 0:
        ugr.append(result[i][1:])

eq_ugr = []
for i in range(len(eq_result)):
    if i != 0:
        eq_ugr.append(eq_result[i][1:])


plt.title('sht #' + str(Shotn) + ', time = ' + str(time) + ' ms')
plt.xlim(0, 1)
plt.ylim(-0.7, 0.7)
plt.grid()
#CP = plt.contour(rgr, zgr, ugr, levels=50)
#plt.clabel(CP, CP.levels)
plt.plot(rb, zb, 'ro-')
for dot in dots.keys():
    plt.scatter(dots[dot][0], dots[dot][1], marker='x')
plt.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
plt.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')


r_lim = [float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])]
z_lim = [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])]
plt.plot(r_lim, z_lim, 'mo-')

strike_point = {'inner': [], 'outer': []}
if min(zb) < dots['x-dot'][1]:
    for j, ind in enumerate([len(rb)-1, 0]):
        for i in range(len(r_lim)):
            if r_lim[i] > rb[ind] and r_lim[i+ 1] < rb[ind]:
                k2, b2 = np.polyfit(r_lim[i:i+2], z_lim[i:i+2], 1)
        print(k2, b2)

        if j:
            ind_first = 0
            ind_last = ind+6
        else:
            ind_first = ind - 5
            ind_last = ind+1
        k1, b1 = np.polyfit(rb[ind_first:ind_last], zb[ind_first:ind_last], 1)

        strike_point[list(strike_point.keys())[j]].extend([(b1-b2)/(k2-k1), (k2*b1-k1*b2)/(k2-k1)])






print(zgr.index(0))

bound = []
for i, e in enumerate(zb):
    if z_need + 1e-3 > e > z_need - 1e-3:
        bound.append(i)
print(bound)

plt.figure()
plt.title('sht #' + str(Shotn) + ', time = ' + str(time) + ' ms')
plt.plot(rgr, ugr[zgr.index(0)], 'x-')
plt.grid()
plt.axvline(rb[bound[0]], color='r')
plt.axvline(rb[bound[1]], color='r')
plt.ylabel('q')
plt.xlabel('R')

plt.figure()
plt.title('sht #' + str(Shotn) + ', time = ' + str(time) + ' ms')
plt.plot(rgr, eq_ugr[zgr.index(0)], 'x-')
plt.grid()
plt.axvline(rb[bound[0]], color='r')
plt.axvline(rb[bound[1]], color='r')
plt.ylabel('psi')
plt.xlabel('R')

"""plt.figure()
plt.plot(zgr, [i[18] for i in ugr])"""
plt.show()
