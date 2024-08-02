import shtRipper
import matplotlib.pyplot as plt
import numpy as np
import json

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_sht_data(Shotn, data_name):
    filename = '//172.16.12.127/Data/sht%i.sht' %Shotn
    res = shtRipper.ripper.read(filename, data_name)
    smooth_res = {}
    print(res.keys())
    for key in res.keys():
        if key != 'nl 42 cm (1.5мм) 64pi':
            baseline = sum(res[key]['y'][:1000]) / len(res[key]['y'][:1000])
            res[key]['y'] = [i - baseline for i in res[key]['y']]
        smooth_res[key] = {}
        smooth_res[key]['data'] = list(smooth(res[key]['y'], 95))
        smooth_res[key]['time'] = res[key]['x']
    return smooth_res

data_name_need = ['Ip внутр.(Пр2ВК) (инт.18)', 'SXR 15 мкм', 'SXR 50 mkm',
                  'nl 42 cm (1.5мм) 64pi', 'D-alfa  хорда R=42 cm',
                  'Emission electrode current', 'Ток пучка новый инжектор']

shot = 43067
data = get_sht_data(shot, data_name_need)
with open('results/%i_sht_data.json' %shot, 'w') as res_file:
    json.dump(data, res_file)

