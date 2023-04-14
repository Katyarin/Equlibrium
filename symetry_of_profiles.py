import json
import matplotlib.pyplot as plt

shotn = 42089

path_res = 'eq_results/%s/' %shotn

with open(path_res + str(shotn) + '_prof_by_psi.json', 'r') as file:
    data_prof = json.load(file)

time = 197.5
P_fig = plt.figure(figsize=(10,5))
p_ax = P_fig.add_subplot(111)
fig = plt.figure()
fig.suptitle(shotn)
Te_ax = fig.add_subplot(121)
ne_ax = fig.add_subplot(122)
now_leg = 'nothing'
count = 0
for t in data_prof.keys():
    if time:
        if float(t) < time - 2 or float(t) > time + 2:
            continue
    print(t)
    #if 165 > float(t) > 153 or 211 > float(t) > 205:
    if 1:
        if float(t) < 160:
            if count == 0:
                fmt_style = 'ro'
            elif count == 1:
                fmt_style = 'r^'
            elif count == 2:
                fmt_style = 'rs'
            else:
                fmt_style = 'r*'
            legend_title_new = 'before L->H (153-160 ms)'
            count+=1
        elif 165 > float(t) > 159:
            if count == 0:
                fmt_style = 'go'
            elif count == 1:
                fmt_style = 'g^'
            elif count == 2:
                fmt_style = 'gs'
            else:
                fmt_style = 'g*'
            legend_title_new = 'after L->H (160-165 ms)'
            count += 1
        else:
            if count == 0:
                fmt_style = 'bo'
            elif count == 1:
                fmt_style = 'b^'
            elif count == 2:
                fmt_style = 'bs'
            else:
                fmt_style = 'b*'
            legend_title_new = 'maximum W (205-211 ms)'
            count += 1
        if now_leg != legend_title_new:
            now_leg, legend_title = legend_title_new, legend_title_new
            count = 0
        else:
            legend_title = None
        ro = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in data_prof[t]['psi']]
        p_ax.scatter(ro, [i * 1000 for i in data_prof[t]['Pe']], label='TS low field side', color='b')
        ro_i = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in
              data_prof[t]['psi_i']]
        p_ax.scatter(ro_i, [i * 1000 for i in data_prof[t]['Pi']], label='CXRS low field side', color='r')

        Te_ax.errorbar(ro, data_prof[t]['Te'], data_prof[t]['Te_err'], fmt=fmt_style, label=legend_title)
        ne_ax.errorbar(ro, data_prof[t]['ne'], data_prof[t]['ne_err'], fmt=fmt_style, label=legend_title)

        p_ax.scatter(ro[-1], 1000 * data_prof[t]['Pe'][-1], marker='o', color='white', label='TS high field side', edgecolors='b')
        p_ax.scatter(ro_i[0], 1000 * data_prof[t]['Pi'][0], marker='o', color='white', label='CXRS high field side', edgecolors='r')
        p_ax.scatter(ro_i[1], 1000 * data_prof[t]['Pi'][1], marker='o', color='white', edgecolors='r')
        ro_prof_e = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in data_prof[t]['prof_e_psi']]
        ro_prof_i = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in data_prof[t]['prof_i_psi']]
        p_ax.plot(ro_prof_e, [i * 1000 for i in data_prof[t]['prof_e_pe']], label='electron', color='b')
        p_ax.plot(ro_prof_i, [i * 1000 for i in data_prof[t]['prof_i_p']], label='ion', color='r')
        Te_ax.errorbar(ro[-1], data_prof[t]['Te'][-1], data_prof[t]['Te_err'][-1], fmt='.', color='white')
        ne_ax.errorbar(ro[-1], data_prof[t]['ne'][-1], data_prof[t]['ne_err'][-1], fmt='.', color='white')

'''for t in data_prof.keys():
    if 'Pi' in list(data_prof[t].keys()):
        ro = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in
                data_prof[t]['psi_P']]
        plt.figure()
        plt.plot(ro, data_prof[t]['P'], label='PET')
        plt.plot(ro, [data_prof[t]['prof_e_pe'][i] + data_prof[t]['prof_i_p'][i] for i in range(len(ro))], label='kinetic')
        plt.title(t)
        plt.grid()
        plt.legend()
plt.show()'''

Te_ax.set_ylabel(r'$T_e$, eV')
ne_ax.set_ylabel(r'$n_e$, $m^{-3}$')
Te_ax.set_xlabel(r'$\rho$')
ne_ax.set_xlabel(r'$\rho$')
p_ax.set_ylabel('P, kPa')
p_ax.set_xlabel(r'$\rho_p$')
Te_ax.legend()
ne_ax.legend()
p_ax.legend()
Te_ax.grid()
ne_ax.grid()
p_ax.grid()
p_ax.set_xlim(0, 1)
p_ax.set_ylim(0, 18)
P_fig.savefig('p_profiles2.png', dpi=600, transparent=True)
plt.show()