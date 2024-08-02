import json
import matplotlib.pyplot as plt

with open('eq_results/statistic.json', 'r') as file:
    stat_data = json.load(file)

need = 'W'

plt.vlines([int(shot) for shot in stat_data[need].keys()],
           [stat_data[need][shot]['av'] - stat_data[need][shot]['range'] for shot in stat_data[need].keys()],
           [stat_data[need][shot]['av'] + stat_data[need][shot]['range'] for shot in stat_data[need].keys()])

for shot in stat_data[need].keys():
    W_max = stat_data[need][shot]['av'] + stat_data[need][shot]['range']
    if W_max > 15000:
        print(W_max, shot, stat_data['shot_without_plasma'][shot])

plt.show()