import pandas
import matplotlib.pyplot as plt
import numpy
import sys
import argparse

def plot_atom(ax, data: pandas.DataFrame, atom: str, ref: pandas.DataFrame):
    subdata = data[data['Atom'] == atom]
    
    lspace = numpy.linspace(-5, 10, 101)
    ax.plot(lspace, lspace, 'k--')
    ax.plot(subdata['Exp'] - ref['Exp'], subdata['SJ'] - ref['SJ'], 'o', label='SJ')
    ax.plot(subdata['Exp'] - ref['Exp'], subdata['SJn'] - ref['SJn'], 'o', label='SJ$^n$')
    
    err_SJ = subdata['Exp'] - ref['Exp'] - (subdata['SJ'] - ref['SJ'])
    err_SJn = subdata['Exp'] - ref['Exp'] - (subdata['SJn'] - ref['SJn'])
    
    ax.text(-4, 8, '{} 1s (N={})'.format(atom, len(err_SJ)), fontsize=12)
    ax.text(-2, 6, '(SJ) {:.2f} $\\pm$ {:.2f}'.format(numpy.mean(err_SJ), numpy.std(err_SJ)), color='tab:blue')
    ax.text(1, -4, '(SJ$^n$) {:.2f} $\\pm$ {:.2f}'.format(numpy.mean(err_SJn), numpy.std(err_SJn)), color='tab:orange')
    
    ax.set_xlim(-5, 10)
        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_C185.csv')
parser.add_argument('-o', '--output', default='Data_XPS_C185.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)

REFS = {  # references
    'B': data.iloc[0],
    'C': data.iloc[35],
    'F': data.iloc[127],
    'N': data.iloc[142],
    'O': data.iloc[162]
}

# print(REFS)

figure = plt.figure(figsize=(7, 10))
(ax1, ax2), (ax3, ax4), (ax5, ax6) = figure.subplots(3, 2, sharey=True)
figure.delaxes(ax6)

ax1.set_ylim(-5, 10)

plot_atom(ax1, data, 'C', REFS['C'])
plot_atom(ax2, data, 'N', REFS['N'])
plot_atom(ax3, data, 'O', REFS['O'])
plot_atom(ax4, data, 'B', REFS['B'])
plot_atom(ax5, data, 'F', REFS['F'])

[ax.set_ylabel('Computed $\\Delta$BE (eV)') for ax in [ax1, ax3, ax5]]
[ax.set_xlabel('Experimental $\\Delta$BE (eV)') for ax in [ax4, ax5]]

plt.tight_layout()
figure.savefig(args.output)
