import json
import matplotlib.pyplot as plt

shotn = 42705

path_res = 'eq_results/%s/' %shotn

with open(path_res + str(shotn) + '_prof_by_psi_0.json', 'r') as file:
    data_prof = json.load(file)

P_fig = plt.figure()
p_ax = P_fig.add_subplot(111)
fig = plt.figure()
fig.suptitle(shotn)
Te_ax = fig.add_subplot(121)
ne_ax = fig.add_subplot(122)
now_leg = 'nothing'
count = 0
for t in data_prof.keys():
    if 165 > float(t) > 153 or 211 > float(t) > 205:
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
        ro = [(i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']) for i in data_prof[t]['psi']]
        p_ax.scatter(ro, data_prof[t]['Pe'])

        Te_ax.errorbar(ro, data_prof[t]['Te'], data_prof[t]['Te_err'], fmt=fmt_style, label=legend_title)
        ne_ax.errorbar(ro, data_prof[t]['ne'], data_prof[t]['ne_err'], fmt=fmt_style, label=legend_title)

        p_ax.plot(ro[-1], data_prof[t]['Pe'][-1], '.', color='white')
        Te_ax.errorbar(ro[-1], data_prof[t]['Te'][-1], data_prof[t]['Te_err'][-1], fmt='.', color='white')
        ne_ax.errorbar(ro[-1], data_prof[t]['ne'][-1], data_prof[t]['ne_err'][-1], fmt='.', color='white')
Te_ax.set_ylabel(r'$T_e$, eV')
ne_ax.set_ylabel(r'$n_e$, $m^{-3}$')
Te_ax.set_xlabel(r'$\rho$')
ne_ax.set_xlabel(r'$\rho$')
Te_ax.legend()
ne_ax.legend()
Te_ax.grid()
ne_ax.grid()
plt.show()