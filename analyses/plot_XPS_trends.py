import pandas
import matplotlib.pyplot as plt
import numpy
import argparse
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

def plot_trends(ax, data: pandas.DataFrame, refs: pandas.DataFrame, labels: list):
    
    for i, label in enumerate(labels):
        subdata_lab = data[data['Type'] == label]
        mi, ma = subdata_lab['SJn'].min(), subdata_lab['SJn'].max()
        
        ax.barh(i, width=ma - mi, height=0.5, left=mi - refs[label[0]]['SJn'], color='black')
    

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_C185.csv')
parser.add_argument('-o', '--output', default='Data_XPS_trends.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)

REFS = {  # references
    'B': data[(data['Molecule'] == 'molecule_001') & (data['Atom' ] == 'B')].iloc[0],
    'C': data[(data['Molecule'] == 'molecule_011') & (data['Atom' ] == 'C')].iloc[0],
    'F': data[(data['Molecule'] == 'molecule_055') & (data['Atom' ] == 'F')].iloc[0],
    'N': data[(data['Molecule'] == 'molecule_018') & (data['Atom' ] == 'N')].iloc[0],
    'O': data[(data['Molecule'] == 'molecule_036') & (data['Atom' ] == 'O')].iloc[0]
}

figure = plt.figure(figsize=(7, 5))
ax1 = figure.subplots(1, 1)

plot_trends(ax1, data, REFS, [
    'C3-C', 'C2-C', 'C1-C', 
    'C3-N', 'C2-N', 'C1-N',
    'C3-O', 'C2-O',
    'C3-F', 'C2-F',
    'N3', 'N2', 'N1',
    'O3', 'O2',
    'F',
    'B',
]
)
ax1.set_xlabel('$\\Delta$BE (eV)')

true_labels =  [
    'C$_{sp^3}$-C', 'C$_{sp^2}$-C', 'C$_{sp}$-C',
    'C$_{sp^3}$-N', 'C$_{sp^2}$-N', 'C$_{sp}$-N',
    'C$_{sp^3}$-O', 'C$_{sp^2}$-O', 
    'C$_{sp^3}$-F', 'C$_{sp^2}$-F',
    'N$_{sp^3}$', 'N$_{sp^2}$', 'N$_{sp}$',
    'O$_{sp^3}$', 'O$_{sp^2}$',
    'F',
    'B',
]
ax1.set_yticks(numpy.arange(len(true_labels)), labels=true_labels)
ax1.invert_yaxis()

ax1.set_xlim(-4, 12)
mi, ma = ax1.get_xlim()
for sep in [2, 5, 7, 9, 12, 14]:
    ax1.plot([mi, ma], [sep  + .5, sep + .5], '--', color='grey')
    
ax1.invert_xaxis()

plt.tight_layout()
figure.savefig(args.output)
