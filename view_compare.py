import matplotlib.pyplot as plt
import json
i_need = ['Ipf1', 'Ipf3', 'Ics', 'Icc', 'Ihfc', 'Ivfc', 'Ipl']


def MAPE(y, aver_y):
    sum_err = 0
    if abs(max(y)) > abs(min(y)):
        y_max = abs(max(y))
    else:
        y_max = abs(min(y))
    for i in range(len(y)):
        if abs(y[i]) >  y_max * 0.04:
            sum_err += abs((abs(y[i]) - abs(aver_y[i])) / abs(y[i]))
    mape = 100 * sum_err / len(y)
    return mape


with open('shot_compare_100.json', 'r') as file:
    all_diff = json.load(file)

Ipf1 = {'niifa': {'small': [], 'big': []}, 'globus': {'small': [], 'big': []}}
for shot in range(11):
    small1 = []
    big1 = []
    small2 = []
    big2 = []
    for i in range(len(all_diff['globus']['Ipf1'][shot])):
        if abs(all_diff['niifa']['Ipf1'][shot][i]) < 2.4:
            small1.append(all_diff['globus']['Ipf1'][shot][i])
            small2.append(all_diff['niifa']['Ipf1'][shot][i])
        else:
            big1.append(all_diff['globus']['Ipf1'][shot][i])
            big2.append(all_diff['niifa']['Ipf1'][shot][i])
    Ipf1['niifa']['small'].append(small2)
    Ipf1['niifa']['big'].append(big2)
    Ipf1['globus']['small'].append(small1)
    Ipf1['globus']['big'].append(big1)


for I in i_need:
    print('____________')
    print(I)
    plt.figure()
    plt.title(I)
    plt.xlabel('niifa, kA')
    plt.ylabel('globus, kA')
    for shot in range(len(all_diff['shots'])):
        mid = int(len(all_diff['globus'][I][shot]) / 2)
        if max(all_diff['niifa']['Ics'][shot]) < 10:
            continue
        if (all_diff['niifa'][I][shot][mid] * all_diff['globus'][I][shot][mid]) / abs(all_diff['niifa'][I][shot][mid] * all_diff['globus'][I][shot][mid]) < 0:
            factor = -1
        else:
            factor = 1
        #plt.plot(all_diff['niifa'][shot][I], [i * factor for i in ['niifa'][I][shot]], 'r--')
        if all_diff['shots'][shot] != 38715:
            plt.scatter(all_diff['niifa'][I][shot], all_diff['globus'][I][shot], s=10, cmap='jet', label=all_diff['shots'][shot])
            #print('--', all_diff['shots'][shot])
            #print(MAPE(all_diff['globus'][I][shot], all_diff['niifa'][I][shot]), '%')
            '''if I == 'Ipf1':
                print('for small A')
                print(MAPE(Ipf1['globus']['small'][shot], Ipf1['niifa']['small'][shot]), '%')
                print('for big A')
                print(MAPE(Ipf1['globus']['big'][shot], Ipf1['niifa']['big'][shot]),  '%')'''
    plt.grid()
    #plt.legend(title='shot #')
    plt.savefig('plots/' + I + '.png', dpi=600)

plt.show()