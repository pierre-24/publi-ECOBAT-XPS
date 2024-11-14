import pandas
import matplotlib.pyplot as plt
import numpy
import argparse


def prepare_data(data: pandas.DataFrame, data_ref: pandas.DataFrame):
    # find corresponding computed data
    computed = []
    for line in data_ref.itertuples():
        subdata = data[data['System'] == line.Molecule]
        subdata = subdata[subdata['Atom_indices'].str.contains(line.Atom)]
        computed.append(subdata.iloc[0]['Value'])
        
    final_data = data_ref.copy()
    final_data.insert(3, 'Computed', computed)
    
    REFS = {  # references
        'B': final_data[(final_data['Molecule'] == 'molecule_001') & (final_data['Atom' ] == 'B_001')].iloc[0],
        'C': final_data[(final_data['Molecule'] == 'molecule_011') & (final_data['Atom' ] == 'C_001')].iloc[0],
        'F': final_data[(final_data['Molecule'] == 'molecule_055') & (final_data['Atom' ] == 'F_001')].iloc[0],
        'N': final_data[(final_data['Molecule'] == 'molecule_018') & (final_data['Atom' ] == 'N_001')].iloc[0],
        'O': final_data[(final_data['Molecule'] == 'molecule_036') & (final_data['Atom' ] == 'O_001')].iloc[0]
    }
    
    delta_exp = []
    delta_computed = []
    
    for line in final_data.itertuples():
        atom = line.Atom[0]
        delta_exp.append(line.Exp - REFS[atom]['Exp'])
        delta_computed.append(line.Computed - REFS[atom]['Computed'])
    
    final_data.insert(4, 'Delta_exp', delta_exp)
    final_data.insert(5, 'Delta_computed', delta_computed)
    
    return final_data
        

def plot_atom(ax, data: pandas.DataFrame, atom: str, color: str, position: tuple, label: str):
    subdata = data[data['Atom'].str.contains(atom)]
    
    ax.plot(subdata['Delta_exp'], subdata['Delta_computed'], 'o', color=color)
    error = subdata['Delta_exp'] - subdata['Delta_computed']
    ax.text(*position, '{:.2f} $\\pm$ {:.2f} ({})'.format(numpy.mean(error), numpy.std(error), label), color=color, rotation=45)
    
    if color == 'tab:blue':
        ax.text(-4, 8, '{} 1s (N={})'.format(atom, len(error)), fontsize=12)
        
parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-r', '--ref', default='../data/Data_XPS_C185_exp.csv')
parser.add_argument('-o', '--output', default='Data_XPS_C185.pdf')

args = parser.parse_args()

if len(args.inputs) != len(args.names):
    raise Exception('len({}) and len({}) do not match'.format(repr(args.inputs), repr(args.names)))

data_ref = pandas.read_csv(args.ref)

data = []
for inp in args.inputs:
    data.append(prepare_data(pandas.read_csv(inp), data_ref))

figure = plt.figure(figsize=(6, 9))
axes = figure.subplots(3, 2, sharey=True)
figure.delaxes(axes[2, 1])

lspace = numpy.linspace(-5, 10, 101)

[ax.plot(lspace, lspace, 'k--') for ax in axes.flatten()]
[ax.set_xlim(-5, 10) for ax in axes.flatten()]
[ax.set_ylim(-5, 10) for ax in axes.flatten()]
[ax.set_ylabel('Computed $\\Delta$BE (eV)') for ax in [axes[0, 0], axes[1, 0], axes[2, 0]]]
[ax.set_xlabel('Experimental $\\Delta$BE (eV)') for ax in [axes[1, 1], axes[2, 0]]]

COLORS = ['tab:blue', 'tab:green']
POSITIONS = [(-4.5, -1), (-1, -4.5)]

for i, subdata in enumerate(data):
    plot_atom(axes[0, 0], subdata, 'C', COLORS[i], POSITIONS[i], args.names[i])
    plot_atom(axes[0, 1], subdata, 'N', COLORS[i], POSITIONS[i], args.names[i])
    plot_atom(axes[1, 0], subdata, 'O', COLORS[i], POSITIONS[i], args.names[i])
    plot_atom(axes[1, 1], subdata, 'B', COLORS[i], POSITIONS[i], args.names[i])
    plot_atom(axes[2, 0], subdata, 'F', COLORS[i], POSITIONS[i], args.names[i])

plt.tight_layout()
figure.savefig(args.output)

