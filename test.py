import Find_boundary
import Globus_PET
import json

result = []
result2 = []
rb = []
zb = []
#PATH = 'c:/work/equilibrium/Globus_PET/'
Path_res = 'results/'

Shotn = 42123
time = 192.5 #ms

with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_eq_res.txt', 'r') as file:
    for line in file:
        result.append([float(i) for i in line.split()])

with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_ipr_res.txt', 'r') as file3:
    for line in file3:
        result2.append([float(i) for i in line.split()])

with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_bound.txt', 'r') as file2:
    for line in file2:
        data = line.split()
        rb.append(float(data[0]))
        zb.append((float(data[1])))

dots = {}
with open(Path_res + str(Shotn) + '_' + str(round(time/1e3, 3)) + '_dots.txt', 'r') as file3:
    for line in file3:
        data = line.split()
        if len(data) < 4:
            dots[data[0]] = [float(data[1]), float(data[2])]
        else:
            dots[data[0] + '_' + data[1]] = [float(data[2]), float(data[3])]
print(dots)
rm = dots['magnetic_axis'][0]

zgr = result[0]
rgr = [result[i][0] for i in range(1, len(result))]
ugr = []
ipr = []
for i in range(len(result)):
    if i != 0:
        ugr.append(result[i][1:])

for i in range(len(result2)):
    if i != 0:
        ipr.append(result2[i][1:])


We_el, error, W_ion, ne_av = Globus_PET.We(Shotn, time/1e3, rgr, zgr, ugr, ipr, rm, Wi=True)
print('result:')
print(We_el, W_ion)
print(We_el+W_ion)
print(ne_av)