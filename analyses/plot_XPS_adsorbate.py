import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import math

from XPS.commons import create_spectrum

FAC = 2.2

ADS = [
    ('C2H6', 'tab:blue'), 
    ('C3H8', 'tab:red'), 
    ('C2H4', 'tab:green'), 
    ('C3H6', 'tab:pink'), 
    ('CO2', 'tab:orange'), 
    ('CH2O', 'tab:orange'), 
    ('C2H4O', 'tab:orange'), 
    ('THF', 'tab:orange'),
    ('CH4O', 'tab:orange'),
    ('THFn', 'tab:orange'),
]

def plot(ax, data: pandas.DataFrame, atom: str, ref: float, xrange: tuple):
    x = numpy.linspace(*xrange, 501)
    
    mins, maxes = [], []
    
    for i, (adsorbate, color) in enumerate(ADS):
        subdata = data[(data['System'] == '{}@Ca'.format(adsorbate)) & (data['Atom'] == atom)]
        
        mins.append(subdata['BE'].min() - ref)
        maxes.append(subdata['BE'].max() - ref)
        
        spectrum = create_spectrum(subdata, x, ref)
        ax.plot(x, FAC * i + spectrum, '-', color='tab:blue')
            
    
    mi, ma = min(filter(lambda x: not numpy.isnan(x), mins)) - 1, max(filter(lambda x: not numpy.isnan(x), maxes)) + 1
    ax.set_xlim(mi, ma)
    
    ax.text(ma - .1, 9 * FAC + 1.7, '{} 1s'.format(atom), fontsize=18)
    
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slab_adsorbate_ref_SJ.csv')
parser.add_argument('-t', '--type', default='SJ')

parser.add_argument('-o', '--output', default='Data_XPS_adsorbate.svg')

args = parser.parse_args()

data = pandas.read_csv(args.input)

if args.type == 'SJn':
    REF_Ca = 426.56
    REF_C = 297.04
    REF_O = 548.79 # in SJn
else:
    REF_Ca = 426.53
    REF_C = 295.31
    REF_O = 546.08 # SJ

figure = plt.figure(figsize=(6, 10))
ax1, ax2 = figure.subplots(1, 2, sharey=True)

plot(ax1, data, 'C', REF_C, (-10, 10))
plot(ax2, data, 'O', REF_O, (-10, 10))

[ax.invert_xaxis() for ax in [ax1, ax2]]
[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax1, ax2]]
[ax.yaxis.set_major_formatter('') for ax in [ax1, ax2]]
[ax.xaxis.set_major_formatter('{x:.1f}') for ax in [ax1, ax2]]

plt.tight_layout()
figure.savefig(args.output)
