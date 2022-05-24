import pandas as pd
import json
import matplotlib.pyplot as plt


with open('data/segment_HA.SLAC.json', 'r') as f:
    res_MEME = json.loads(f.read())
print(res_MEME.keys())


for k, v in res_MEME['analysis'].items():
    print(k.capitalize(), v, sep=': ')

print(res_MEME['MLE'].keys())
print(res_MEME['branch attributes'].keys())
# print(res_MEME['data partitions'])

# for k, v in res_MEME['fits'].items():
#     for k_, v_ in v.items():
#         print(k, k_, v_)

print(res_MEME['runtime'])
print(res_MEME['sample-2.5']['0'].keys())
print(res_MEME['sample-97.5']['0'].keys())
print(res_MEME['sample-median']['0'].keys())
# print(res_MEME['tested']['0'].keys())

# for key_1 in res_MEME.keys():
#     if key_1 not in ['tested', 'sample-median', 'sample-2.5', 'sample-97.5', 'input', 'fits']:
#         pass
        # print(key_1, res_MEME[key_1])
