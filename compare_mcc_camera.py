import json
from scipy import interpolate as inter
import matplotlib.pyplot as plt

with open('camera_data.json', 'r') as file:
    camera = json.load(file)

mcc_data = {}
shot_list = []
p = 0
with open('mcc_str_point.txt', 'r') as file2:
    for line in file2:
        line_data = line.split()
        if p == 0:
            for element in line_data:
                mcc_data[element] = []
                shot_list.append(element)
        else:
            for i, el in enumerate(line_data):
                mcc_data[shot_list[i]].append(float(el))
        p+=1

print(camera.keys())
print(mcc_data.keys())
comparison = {}

for shot in camera.keys():
    comparison[shot] = {'camera': [], 'mcc': []}
    min_time = min(camera[shot]['t'])
    max_time = max(camera[shot]['t'])
    print(min_time, max_time)
    camera_int = inter.interp1d(camera[shot]['t'], camera[shot]['data'], kind='linear')
    for el, time in enumerate(mcc_data['Time']):
        print(time, el)
        if max_time > time > min_time:
            comparison[shot]['camera'].append(camera_int(time))
            comparison[shot]['mcc'].append(mcc_data[shot][el] * 10)

plt.figure()

MAPE_by_shot = []

for shot in comparison.keys():
    mape = 0
    for i in range(len(comparison[shot]['camera'])):
        mape += abs((comparison[shot]['mcc'][i] - comparison[shot]['camera'][i]) / comparison[shot]['camera'][i])
    MAPE_by_shot.append(mape/ len(comparison[shot]['camera']))

print(shot_list)
print(MAPE_by_shot)

plt.plot(shot_list[1:], MAPE_by_shot, 'o')
plt.grid()

plt.figure()

for shot in comparison.keys():
    plt.scatter(comparison[shot]['camera'], comparison[shot]['mcc'], label=shot)
plt.plot([i for i in range(200, 400)], [i for i in range(200, 400)])
plt.xlim(200, 400)
plt.ylim(200, 400)
plt.legend()
plt.show()
