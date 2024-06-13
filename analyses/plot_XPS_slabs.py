import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse

def gaussian(x, mu, sigma=0.2):
    return 1/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-.5*((x-mu)/sigma)**2)

def plot_slab(ax, data: pandas.DataFrame, system: str, atom: str, ref: float, xrange: tuple):
    lspace = numpy.linspace(*xrange, 201)
    
    for i, subsystem in enumerate([3, 4, 5, 6, 7, 8]):
        if subsystem == 8 and 'CaH2' in system:
            continue
        
        subdata = data[data['System'] == '{}/{}'.format(system, subsystem)]
        subdata = subdata[subdata['Atom'] == atom]
        
        yspace = numpy.zeros(lspace.shape)
        N = 0
        for BE, n in zip(subdata['BE'], subdata['N']):
            yspace += gaussian(lspace, BE - ref) * n
            N += n
        
        ax.plot(lspace, i + yspace / N, '-')
        ax.text(xrange[1] - .25, i + .25, subsystem * (4 if system == 'CaH2' else 2))
    
    ax.text(xrange[0], 6.5, '{} / {}{}s'.format(system.replace('H2', 'H$_2$'), atom, 2 if atom == 'Ca' else 1), fontsize=12)

        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_slabs.csv')
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)
# REF_CA = 426.53
# REF_O = 546.08 # SJ
REF_CA = 426.56
REF_O = 548.79 # in SJn

figure = plt.figure(figsize=(10, 7))
(ax1, ax2), (ax3, ax4) = figure.subplots(2, 2, sharey=True)

plot_slab(ax1, data, 'Ca', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax2, data, 'CaH2', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax3, data, 'CaO', 'Ca', REF_CA, (-1.5, 1.5))
plot_slab(ax4, data, 'CaO', 'O', REF_O, (-9, -6))

[ax.set_xlabel('$\\Delta$BE (eV)') for ax in [ax3, ax4]]

plt.tight_layout()
figure.savefig(args.output)
