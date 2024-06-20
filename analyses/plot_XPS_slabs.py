import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse

from XPS.commons import create_spectrum

def plot_slab(ax, data: pandas.DataFrame, system: str, atom: str, ref: float, xrange: tuple):
    x = numpy.linspace(*xrange, 201)
    
    for i, subsystem in enumerate([3, 4, 5, 6, 7, 8]):
        if subsystem == 8 and 'CaH2' in system:
            continue
        
        subdata = data[(data['System'] == '{}/{}'.format(system, subsystem)) & (data['Atom'] == atom)]
        
        ax.plot(x, i + create_spectrum(subdata, x, ref), '-')
        ax.text(xrange[1] - .25, i + .25, subsystem * (4 if system == 'CaH2' else 2))
    
    ax.text(xrange[1], 6.5, '{} / {}{}s'.format(system.replace('H2', 'H$_2$'), atom, 2 if atom == 'Ca' else 1), fontsize=12)

        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slab_SJ.csv')
parser.add_argument('-t', '--type', default='SJ')
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)

if args.type == 'SJn':
    REF_CA = 426.56
    REF_O = 548.79 # in SJn
else:
    REF_CA = 426.53
    REF_O = 546.08 # SJ

figure = plt.figure(figsize=(10, 7))
(ax1, ax2), (ax3, ax4) = figure.subplots(2, 2, sharey=True)

plot_slab(ax1, data, 'Ca', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax2, data, 'CaH2', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax3, data, 'CaO', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax4, data, 'CaO', 'O', REF_O, (-9, -6) if args.type == 'SJn' else (-7, -4))

[ax.invert_xaxis() for ax in [ax1, ax2, ax3, ax4]]
[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax3, ax4]]
[ax.yaxis.set_major_formatter('') for ax in [ax1, ax3]]

plt.tight_layout()
figure.savefig(args.output)
