import Find_boundary
import numpy as np
import Globus_PET
import subprocess
import matplotlib.pyplot as plt
import time as t
from matplotlib.backends.backend_pdf import PdfPages

pi = 3.14159265359
mu0 = 4*pi*1e-7
PATH = 'c:/work/equilibrium/Globus_PET/'

Shotn = 38800
time_list = [0.181]

betta_p_list = []
betta_I_list = []
li_find_list = []

pdf_file = PdfPages(str(Shotn) + 'li_bp_plots.pdf')

for ind in range(len(time_list)):
    print('______________________________________________________________________________________')
    print(ind)
    time = time_list[ind]
    f = '//172.16.12.127/Pub/!!!CURRENT_COIL_METHOD/old_mcc/mcc_%d.json' % Shotn
    time, Rc, Rcm = Find_boundary.bound(f, time)
    print('new time: ', time)
    I_coil = Globus_PET.get_coils(Shotn, time)
    Globus_PET.COIL_upd(Shotn, time, I_coil)

    print(I_coil['Ipl'], Rc, Rcm)

    for betta_I in range(1, 9):
        betta_I = betta_I/10
        print('betta_I: ', betta_I)
        betta_I_list.append(betta_I)
        li = 1
        #alf11, alf22 = 1.5, 1
        li_f, li_f2, res = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [50, 190, 10], pdf_file, True, 'all')
        #li_prob = int(li_f2 * 100)
        #li_acc, li_acc2, res_acc = Globus_PET.find_par2('li', Shotn, time, I_coil, betta_I, li, [li_prob - 7, li_prob + 7, 1], pdf_file, True, 'all')
        #plt.show()
        #print('real li: ', res['li'])
        #print('real bp: ', res['bp'])
        print(li_f)
        betta_p_list.append(res['bp'][res['li'].index(li_f)])
        li_find_list.append(li_f)

pdf_file.close()