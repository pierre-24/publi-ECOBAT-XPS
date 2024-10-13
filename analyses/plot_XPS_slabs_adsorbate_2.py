import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import math
from scipy.signal import argrelextrema

from XPS.commons import create_spectrum_BE

FAC = 3.75
SLBS = [('Ca', 'tab:blue'), ('CaO', 'tab:red'), ('CaO_OH2', 'tab:green'), ('CaH2', 'tab:pink')]


def prepare_data(data: pandas.DataFrame, no_C: bool = True):
    REFS = {  # references
        'Ca': data[(data['System'] == 'ref_Ca') | ((data['System'] == 'CaO_slab/8') & (data['Atom_indices'].str.contains('Ca_025')))].iloc[0],
        'O': data[(data['System'] == 'ref_O')].iloc[0],
        'C': -1 if no_C else data[(data['System'] == 'ref_C')].iloc[0]
    }
    
    delta_computed = []
    data_out = data[~data['System'].isin(['ref_Ca', 'ref_O'])]
    
    for line in data_out.itertuples():
        atom = line.Atom
        delta_computed.append(line.Value - REFS[atom]['Value'])
    
    data_out.insert(5, 'Delta_computed', delta_computed)
    
    return data_out

def annotate(ax, data: pandas.DataFrame, annotations: list):
    pass
    
def get_annotations(inp: str):
    annotations = {}
    for annotation in inp.split(','):
        df,label = annotation.split('=')
        atom, system = df.split('@')
        symbol, _ = atom.split('_')
        
        if symbol not in annotations:
            annotations[symbol] = {}
        
        if system not in annotations[symbol]:
            annotations[symbol][system] = []
            
        annotations[symbol][system].append((atom, label))
    
    return annotations
        

def plot_atom(ax, data: pandas.DataFrame, data_slabs: pandas.DataFrame, adsorbate: str, atom: str, xrange: tuple, annotations: dict):
    x = numpy.linspace(*xrange, 501)
    
    mins, maxes = [], []
    
    for i, (slab, color) in enumerate(SLBS):
        slab_system = '{}_slab/3'.format(slab)
        slab_ads_system = '{}_{}'.format(slab, adsorbate)
        
        subdata = data[(data['System'] == slab_ads_system) & (data['Atom'] == atom)]
        subdata_slabs = data_slabs[(data_slabs['System'] == slab_system) & (data_slabs['Atom'] == atom)]
        
        ax.plot(xrange, [i * FAC, i * FAC], '-', color='grey')
        
        y = create_spectrum_BE(subdata, x)
        ax.plot(x, i * FAC + y, color=color)
        
        if atom != 'C':
            y_slab = create_spectrum_BE(subdata_slabs, x)
            ax.plot(x, i * FAC - y_slab, '--', color=color)
            ax.plot(x, i * FAC - y_slab + y, ':', color=color)
        
        mins.extend([subdata['Delta_computed'].min(), subdata_slabs['Delta_computed'].min()])
        maxes.extend([subdata['Delta_computed'].max(), subdata_slabs['Delta_computed'].max()])
        
        if slab_ads_system in annotations:
            for a_atom, a_label in annotations[slab_ads_system]:
                d = subdata[subdata['Atom_indices'].str.contains(a_atom)]
                if len(d) > 0:
                    v = d.iloc[0]['Delta_computed']
                    iloc = int((v - xrange[0]) / (xrange[1] - xrange[0]) * 501)
                    ax.text(v, i * FAC + y[iloc] + .1, a_label, va='bottom', ha='center', color=color)
                
        if slab_system in annotations:
            for a_atom, a_label in annotations[slab_system]:
                d = subdata_slabs[subdata_slabs['Atom_indices'].str.contains(a_atom)]
                if len(d) > 0:
                    v = d.iloc[0]['Delta_computed']
                    iloc = int((v - xrange[0]) / (xrange[1] - xrange[0]) * 501)
                    ax.text(v, i * FAC - y_slab[iloc] - .1, a_label, va='top', ha='center', color=color)
                
    
    mi, ma = min(filter(lambda x: not numpy.isnan(x), mins)) - 1, max(filter(lambda x: not numpy.isnan(x), maxes)) + 1
    ax.set_xlim(mi, ma)
    
    ax.text(.05, .95, '{} {}s'.format(atom, 2 if atom == 'Ca' else 1), fontsize=18, transform=ax.transAxes)
    
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/results_slab_adsorbate/Data_XPS_slabs_adsorbate_SJ_Evac.csv')
parser.add_argument('-is', '--input-slabs', default='../data/results_slabs/Data_XPS_slabs_SJ_Evac.csv')
parser.add_argument('-a', '--adsorbate', default='C2H4')
parser.add_argument('-n', '--name', default='C2H4')
parser.add_argument('-x', '--annotate', type=get_annotations)
parser.add_argument('-o', '--output', default='Data_XPS_slab_adsorbate.pdf')

args = parser.parse_args()

data = prepare_data(pandas.read_csv(args.input), no_C=False)
data_slabs = prepare_data(pandas.read_csv(args.input_slabs))

figure = plt.figure(figsize=(12, 10))
axes = figure.subplots(1, 3, sharey=True)
axes[0].set_ylim(-3, FAC * 4 + .5)

annotations = args.annotate if args.annotate is not None else {}

plot_atom(axes[0], data, data_slabs, args.adsorbate, 'Ca', (-5, 5), args.annotate['Ca'] if 'Ca' in annotations else {})
plot_atom(axes[1], data, data_slabs, args.adsorbate, 'O', (-10, 10), args.annotate['O'] if 'O' in annotations else {})
plot_atom(axes[2], data, data_slabs, args.adsorbate, 'C', (-5, 5), args.annotate['C'] if 'C' in annotations else {})

[ax.invert_xaxis() for ax in axes]
[ax.set_xlabel('$\\Delta$BE (eV)') for ax in axes]
[ax.yaxis.set_major_formatter('') for ax in axes]
[ax.xaxis.set_major_formatter('{x:.1f}') for ax in axes]
[ax.tick_params('y', left=False, labelleft=False) for ax in axes]

axes[1].text(.5, 1.05, args.name, fontsize=18, horizontalalignment='center', transform=axes[1].transAxes)

plt.tight_layout()
figure.savefig(args.output)
