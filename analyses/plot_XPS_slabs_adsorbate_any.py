import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse
import math
from scipy.signal import argrelextrema

from XPS.commons import create_spectrum

FAC = 0.7

def plot_atom(ax, data: pandas.DataFrame, systems: list, atom: str, ref: float, xrange: tuple):
    x = numpy.linspace(*xrange, 501)
    
    for i, system in enumerate(systems):
        subdata = data[(data['System'] == system) & (data['Atom'] == atom)]
        
        # make spectrum
        spectrum = create_spectrum(subdata, x, ref, sigma=0.6)
        ax.plot(x, FAC * i + spectrum, '-')
        ax.plot(xrange, [FAC * i, FAC * i], color='grey')
        
    ax.text(xrange[1] - .1, len(systems) * FAC - .2, '{} {}s'.format(atom, 2 if atom == 'Ca' else 1), fontsize=18)
    ax.set_xlim(*xrange)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slab_adsorbate_SJ.csv')
parser.add_argument('-t', '--type', default='SJ')
parser.add_argument('-s', '--systems', default='THF@Ca THFn@Ca THF@CaO THFn@CaO THF@CaO_OH2 THFn@CaO_OH2', type=lambda x: x.split(' '))
parser.add_argument('-n', '--names', default='THF@Ca HO(C$_4$H$_8$O)$_2$H@Ca THF@CaO HO(C$_4$H$_8$O)$_2$H@CaO THF@CaO·H$_2$O HO(C$_4$H$_8$O)$_2$H@CaO·H$_2$O', type=lambda x: x.split(' '))

parser.add_argument('-o', '--output', default='Data_XPS_slab_adsorbate_any.pdf')

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

figure = plt.figure(figsize=(8, 6))
ax2, ax3 = figure.subplots(1, 2, sharey=True)

ax2.set_ylim(0, len(args.systems) * FAC)

plot_atom(ax2, data, args.systems, 'C', REF_C, (-4, 6) if args.type == 'SJ' else (-6,4))
plot_atom(ax3, data, args.systems, 'O', REF_O, (-9, 1) if args.type == 'SJ' else (-11,-1))

for i, name in enumerate(args.names):
    ax2.text(5.75 if args.type == 'SJ' else 3.75, i * FAC + .1, name)

[ax.invert_xaxis() for ax in [ax2, ax3]]
[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax2, ax3]]
[ax.yaxis.set_major_formatter('') for ax in [ax2, ax3]]
[ax.xaxis.set_major_formatter('{x:.1f}') for ax in [ax2, ax3]]

plt.tight_layout()
figure.savefig(args.output)
