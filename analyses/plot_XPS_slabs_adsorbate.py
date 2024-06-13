import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse

from XPS.commons import create_spectrum

FAC = 5

def plot(ax, data: pandas.DataFrame, data_slabs: pandas.DataFrame, data_adsorbates: pandas.DataFrame, adsorbate: str, atom: str, ref: float, xrange: tuple):
    x = numpy.linspace(*xrange, 201)
    subdata_ads = data_adsorbates[(data_adsorbates['System'] == adsorbate) & (data_adsorbates['Atom'] == atom)]
    spectrum_ads = create_spectrum(subdata_ads, x, ref)
    
    for i, (slab, color) in enumerate([('Ca', 'tab:blue'), ('CaO', 'tab:red'), ('CaO_OH2', 'tab:green'), ('CaH2', 'tab:pink')]):
        subdata = data[(data['System'] == '{}@{}'.format(adsorbate, slab)) & (data['Atom'] == atom)]
        subdata_slab = data_slabs[(data_slabs['System'] == '{}/3'.format(slab)) & (data_slabs['Atom'] == atom)]
        
        spectrum = create_spectrum(subdata, x, ref)
        spectrum_slab = create_spectrum(subdata_slab, x, ref)
        
        # TODO: we will probably have to shift the spectrum so that it matches more or less for the innermost atoms (but how to get innermost?)
        # TODO: one molecule on each slab compared to one molecule alone?
        
        ax.plot(x, FAC * i - spectrum_slab - spectrum_ads, '-', color=color)
        ax.plot(x, FAC * i + 1 + spectrum - spectrum_slab - spectrum_ads, '--', color=color)
        ax.plot(x, FAC * i + 2 + spectrum, '-', color=color)
    
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slab_adsorbate.csv')
parser.add_argument('-is', '--input-slabs', default='../data/Data_XPS_slabs.csv')
parser.add_argument('-ia', '--input-adsorbates', default='../data/Data_XPS_adsorbates.csv')
parser.add_argument('-a', '--adsorbate', default='C2H4')

parser.add_argument('-o', '--output', default='Data_XPS_slab_adsorbate.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)
data_slabs = pandas.read_csv(args.input_slabs)
data_adsorbates = pandas.read_csv(args.input_adsorbates)

REF_CA = 426.56
REF_C = 297.04 
REF_O = 548.79 # in SJn

figure = plt.figure(figsize=(10, 7))
ax1, ax2 = figure.subplots(1, 2, sharey=True)

plot(ax1, data, data_slabs, data_adsorbates, args.adsorbate, 'C', REF_C, (-9, 0))
plot(ax2, data, data_slabs, data_adsorbates, args.adsorbate, 'O', REF_O, (-12, -2))

[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax1, ax2]]
[ax.yaxis.set_major_formatter('') for ax in [ax1, ax2]]

plt.tight_layout()
figure.savefig(args.output)
