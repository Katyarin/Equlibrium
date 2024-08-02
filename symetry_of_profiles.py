import json
import matplotlib.pyplot as plt
plt.rc('font', size= 18 )

shotn = 43358

path_res = 'eq_results/%s/' %shotn

with open(path_res + str(shotn) + '_prof_by_psi.json', 'r') as file:
    data_prof = json.load(file)

time = 186
P_fig = plt.figure(figsize=(20,7))
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
        ro_P = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis'])) ** 0.5 for i in
              data_prof[t]['psi_P']]
        p_ax.plot(ro_P, [i * 1e3 for i in data_prof[t]['P']], label='PET all P')
        p_ax.errorbar(ro, [i * 1000 for i in data_prof[t]['Pe']], yerr=[(data_prof[t]['Pe'][i] * (((data_prof[t]['Te_err'][i] / data_prof[t]['Te'][i])**2 + (data_prof[t]['ne_err'][i] / data_prof[t]['ne'][i])**2)**0.5))*1e3 for i in range(len(ro))], label='ТР на стороне слабого магнитного поля', color='b', fmt='o')
        #ro_i = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in
        #      data_prof[t]['psi_i']]
        #p_ax.errorbar(ro_i, [i * 1000 for i in data_prof[t]['Pi']], yerr=[i * 1000 for i in data_prof[t]['Pi_err']], label='CXRS на стороне слабого магнитного поля', color='r', fmt='o')

        Te_ax.errorbar(ro, data_prof[t]['Te'], yerr=data_prof[t]['Te_err'], fmt=fmt_style, label=legend_title)
        ne_ax.errorbar(ro, data_prof[t]['ne'], yerr=data_prof[t]['ne_err'], fmt=fmt_style, label=legend_title)

        p_ax.scatter(ro[-1], 1000 * data_prof[t]['Pe'][-1], marker='o', color='white', label='TP на стороне сильного магнитного поля', edgecolors='b', zorder=10)
        #p_ax.scatter(ro_i[0], 1000 * data_prof[t]['Pi'][0], marker='o', color='white', label='CXRS на стороне сильного магнитного поля', edgecolors='r', zorder=10)
        #p_ax.scatter(ro_i[1], 1000 * data_prof[t]['Pi'][1], marker='o', color='white', edgecolors='r', zorder=10)
        ro_prof_e = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in data_prof[t]['prof_e_psi']]
        #ro_prof_i = [((i - data_prof[t]['psi_axis']) / (data_prof[t]['psi_b'] - data_prof[t]['psi_axis']))**0.5 for i in data_prof[t]['prof_i_psi']]
        p_ax.plot(ro_prof_e, [i * 1000 for i in data_prof[t]['prof_e_pe']], label='давление электронов', color='b')
        #p_ax.plot(ro_prof_i, [i * 1000 for i in data_prof[t]['prof_i_p']], label='давление ионов', color='r')
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

Te_ax.set_ylabel(r'$T_e$, эВ')
ne_ax.set_ylabel(r'$n_e$, $м^{-3}$')
Te_ax.set_xlabel(r'$\rho$')
ne_ax.set_xlabel(r'$\rho$')
p_ax.set_ylabel('P, кПа')
p_ax.set_xlabel(r'$\rho_p$')
Te_ax.legend()
ne_ax.legend()
p_ax.legend()
Te_ax.grid()
ne_ax.grid()
p_ax.grid()
p_ax.set_xlim(0, 1)
p_ax.set_ylim(0, 18)


from PIL import Image
from io import BytesIO

png1 = BytesIO()
P_fig.savefig(png1, format='png')

# (2) load this image into PIL
png2 = Image.open(png1)

# (3) save as TIFF
png2.save('fig3.tif')
png1.close()

P_fig.savefig('p_profiles2.png', dpi=600, transparent=True)
plt.show()