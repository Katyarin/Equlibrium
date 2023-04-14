import Find_boundary
import Globus_PET
import json

tomorow = '230317'
shot_list = [41585, 42368, 42123, 41644, 41645, 41649, 41665, 41666, 41756, 41813, 41820, 41859, 41863,
             41865, 41866, 41867, 42090, 42095, 42097, 42123, 42368, 42414, 42416, 42492,
             42529, 42567, 42705]
W = []
W_approx = []
for Shotn in shot_list:
    data = {}
    path_res = 'eq_results/%s/' % Shotn
    p=0
    list_names = []
    with open(path_res + tomorow + '_res_' + str(Shotn) + '.txt', 'r') as res_file:
        for line in res_file:
            if p==0:
                for i in list(line.split()):
                    data[i] = []
                    list_names.append(i)
                p +=1
            else:
                for i, el in enumerate(list(line.split())):
                    if el:
                        data[list_names[i]].append(float(el))
    W.extend(data['W'])
    W_approx.extend(data['W_approx'])
print(W)
print(W_approx)
MAPE = 0
for i in range(len(W)):
    MAPE += abs(W[i]-W_approx[i])/W[i]
Mape2 =  MAPE/len(W)
print(Mape2)
