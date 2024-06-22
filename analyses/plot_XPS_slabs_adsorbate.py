import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import math

from XPS.commons import create_spectrum

FAC = 3

SLBS = [('Ca', 'tab:blue'), ('CaO', 'tab:red'), ('CaO_OH2', 'tab:green'), ('CaH2', 'tab:pink')]

def plot(ax, data: pandas.DataFrame, data_ref: pandas.DataFrame, data_align: pandas.DataFrame, adsorbate: str, atom: str, ref: float, xrange: tuple):
    x = numpy.linspace(*xrange, 501)
    
    mins, maxes = [], []
    
    for i, (slab, color) in enumerate(SLBS):
        subdata = data[(data['System'] == '{}@{}'.format(adsorbate, slab)) & (data['Atom'] == atom)]
        subdata_ref = data_ref[(data_ref['System'] == '{}@{}'.format(adsorbate, slab)) & (data_ref['Atom'] == atom)]
        
        align = data_align[(data_align['System'] == '{}@{}'.format(adsorbate, slab))].iloc[0]['Ref']
        sa = data[(data['System'] == '{}@{}'.format(adsorbate, slab)) & (data['Atoms'].str.contains(align))].iloc[0]['BE']
        sr = data_ref[(data_ref['System'] == '{}@{}'.format(adsorbate, slab)) & (data_ref['Atoms'].str.contains(align))].iloc[0]['BE']
        iref = sr - sa
        
        mins.extend([subdata['BE'].min() - ref -iref, subdata_ref['BE'].min() - ref])
        maxes.extend([subdata['BE'].max() - ref - iref, subdata_ref['BE'].max() - ref])
        
        spectrum = create_spectrum(subdata, x, ref - iref)
        spectrum_ref = create_spectrum(subdata_ref, x, ref)
        
        ax.plot(x, FAC * i + spectrum_ref, '--', color=color)
        ax.plot(x, FAC * i + spectrum - spectrum_ref, ':', color=color)
        ax.plot(x, FAC * i + spectrum, '-', color=color)
            
    
    mi, ma = min(filter(lambda x: not numpy.isnan(x), mins)) - .5, max(filter(lambda x: not numpy.isnan(x), maxes)) + .5
    ax.set_xlim(mi, ma)
    
    ax.text(ma - .1, 3 * FAC + 2, '{} {}s'.format(atom, 2 if atom == 'Ca' else 1), fontsize=18)
    
    for i, (slab, color) in enumerate(SLBS):
        ax.plot([mi, ma], [FAC * i, FAC * i], color='grey')
        if atom == 'Ca':
            txt = slab
            if txt == 'CaO_OH2':
                txt = 'CaO$\\cdot$H$_2$O'
            elif txt == 'CaH2':
                txt = 'CaH$_2$'
            ax.text(mi + .5, FAC * i + 2, txt, color=color)
    
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slab_adsorbate_SJ.csv')
parser.add_argument('-ir', '--input-ref', default='../data/Data_XPS_slab_adsorbate_ref_SJ.csv')
parser.add_argument('-ia', '--input-align', default='../data/align_SJ.csv')
parser.add_argument('-t', '--type', default='SJ')
parser.add_argument('-a', '--adsorbate', default='C2H4')
parser.add_argument('-n', '--name', default='C2H4')

parser.add_argument('-o', '--output', default='Data_XPS_slab_adsorbate.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)
data_ref = pandas.read_csv(args.input_ref)
data_align = pandas.read_csv(args.input_align)

if args.type == 'SJn':
    REF_Ca = 426.56
    REF_C = 297.04
    REF_O = 548.79 # in SJn
else:
    REF_Ca = 426.53
    REF_C = 295.31
    REF_O = 546.08 # SJ

figure = plt.figure(figsize=(10, 7))
ax1, ax2, ax3 = figure.subplots(1, 3, sharey=True)

plot(ax1, data, data_ref, data_align, args.adsorbate, 'Ca', REF_Ca, (-5, 5))
plot(ax2, data, data_ref, data_align, args.adsorbate, 'C', REF_C, (-10, 10))
plot(ax3, data, data_ref, data_align, args.adsorbate, 'O', REF_O, (-10, 10))

[ax.invert_xaxis() for ax in [ax1, ax2, ax3]]
[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax1, ax2, ax3]]
[ax.yaxis.set_major_formatter('') for ax in [ax1, ax2, ax3]]
[ax.xaxis.set_major_formatter('{x:.1f}') for ax in [ax1, ax2, ax3]]

ax2.text(numpy.mean(ax2.get_xlim()), FAC * 4, args.name, fontsize=18, horizontalalignment='center')

plt.tight_layout()
figure.savefig(args.output)
