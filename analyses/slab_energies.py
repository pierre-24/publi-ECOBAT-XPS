import pandas
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

EV_TO_KJMOL = 96.485 # kJ mol⁻¹

def prepare_data(data: pandas.DataFrame):
    subdata_systems = data[~data['system'].str.contains('@')]
    systems = dict(zip(subdata_systems['system'], subdata_systems['energy']))

    subdata_slabs = data[data['system'].str.contains('@')]
    subdata_slabs.insert(1, 'adsorbate', [x.split('@')[0] for x in subdata_slabs['system']])
    subdata_slabs.insert(1, 'slab', [x.split('@')[1] for x in subdata_slabs['system']])
    
    subdata_slabs.insert(1, 'dE', [(e - 9 * systems[slab] - 2 * systems[adsorbate])* EV_TO_KJMOL for e, slab, adsorbate in zip(subdata_slabs['energy'], subdata_slabs['slab'], subdata_slabs['adsorbate'])])
    
    return subdata_slabs


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_slab_energies.csv')

args = parser.parse_args()

# read data
data = pandas.read_csv(args.input)
data = prepare_data(data)

# group per slab
subdata = {}
for slab in ['Ca', 'CaO', 'CaO_OH2', 'CaH2']:
    subdata[slab] = data[data['slab'] == slab]

table = subdata['Ca'][['adsorbate', 'dE']]\
    .join(subdata['CaO'][['adsorbate', 'dE']].set_index('adsorbate'), on='adsorbate', lsuffix='_Ca', rsuffix='_CaO')\
    .join(subdata['CaO_OH2'][['adsorbate', 'dE']].set_index('adsorbate'), on='adsorbate')\
    .join(subdata['CaH2'][['adsorbate', 'dE']].set_index('adsorbate'), on='adsorbate', lsuffix='_CaO_OH2', rsuffix='_CaH2')

for row in table.iterrows():
    print('{} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\'.format(*row[1]))
