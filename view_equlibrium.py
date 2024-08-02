import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import Find_boundary
import json
result = []
rb = []
zb = []

from PIL import Image
from io import BytesIO

#Path_res = 'triangularity/'
Path_res = 'results/'
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 43922
time = 180 #ms

with open(PATH + 'BLANFW.DAT', 'r') as file3:
    blanfw = []
    for line in file3:
        blanfw.append(line.split())

Rvv = float(blanfw[0][0])
nvv = int(blanfw[1][0])

try:
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/V3_zad7_mcc/mcc_%d.json' % Shotn
    t, Rc, Rcm, Rav, k, conf = Find_boundary.bound(f, time / 1000)
except FileNotFoundError:
    print('not found in new version')
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    t, Rc, Rcm, Rav, k, conf = Find_boundary.bound(f, time / 1000)





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


with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_eq_res.txt', 'r') as file:
    for line in file:
        result.append([float(i) for i in line.split()])

fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax.set_xlim(0, 1)
ax.set_ylim(-0.8, 0.8)
'''for dop in ['_min', '_max']:
    rb = []
    zb = []
    with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_bound' + dop +'.txt', 'r') as file2:
        for line in file2:
            data = line.split()
            rb.append(float(data[0]))
            zb.append((float(data[1])))
    ax.plot(rb, zb, '--', label=('PET' + dop))
rb = []
zb = []'''
with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) + '_bound' +'.txt', 'r') as file2:
    for line in file2:
        data = line.split()
        rb.append(float(data[0]))
        zb.append((float(data[1])))
#ax.plot(rb, zb, label=('PET'))

dots = {}

with open(Path_res + str(Shotn) + '_' + str(round((time/1e3),3)) +'_dots.txt', 'r') as file3:
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

ax.contour(rgr, zgr, ugr, levels=100, alpha=0.5)
ax.plot(rb, zb, 'g', label='PET')

ax.plot([i / 100 for i in mcc_bound['r']], [i / 100 for i in mcc_bound['z']], 'r', label='MMF')

for dot in dots.keys():
    ax.scatter(dots[dot][0], dots[dot][1], marker='x',  zorder=2)
    print(dots[dot][0], dots[dot][1])

ax.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
ax.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')



ax.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
            [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')
#ax.grid()
ax.set_ylabel('z, м')
ax.set_xlabel('r, м')
ax.legend()
ax.set_title('Разряд #%i, время = %i мс' %(Shotn, time) )
#ax.set_facecolor('#DDE9EC')
coils_name = ['pf1', 'pf1', 'pf2', 'pf2', 'pf3', 'pf3'] + ['hfc'] *4 + ['vfc'] * 2 + ['cc'] * 6
print(len(r_coil))
for i in range(len(r_coil)):
    ax.add_patch(Rectangle((r_coil[i] - dr_coil[i]/2, z_coil[i]- dz_coil[i]/2), dr_coil[i], dz_coil[i], ec='blue'))
    if i < 10:
        if i % 2 == 0:
            ax.text(r_coil[i] - dr_coil[i] / 2, z_coil[i] + 0.05, coils_name[i], fontsize=12)
        else:
            ax.text(r_coil[i] - dr_coil[i] / 2, z_coil[i] + 0.05 * -2, coils_name[i], fontsize=12)
    elif i < 18:
        ax.text(r_coil[i] + dr_coil[i] / 2, z_coil[i], coils_name[i], fontsize=12)
#ax.scatter(r_coil, z_coil, zorder=2)

loop_r = []
loop_z = []
with open('loop.txt', 'r') as loop_file:
    for line in loop_file:
        data = line.split()
        loop_r.append(float(data[0]))
        loop_z.append(float(data[1]))

ax.scatter(loop_r, loop_z, color='g')
plt.savefig('plots/%i_%.2f.png' %(Shotn, time), dpi=900, transparent= True)

png1 = BytesIO()
fig.savefig(png1, format='png')

# (2) load this image into PIL
png2 = Image.open(png1)

# (3) save as TIFF
png2.save('fig1.tiff')
png1.close()
plt.show()