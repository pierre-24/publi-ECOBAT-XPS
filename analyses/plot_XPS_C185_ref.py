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
    
    lspace = numpy.linspace(-5, 10, 101)
    ax.plot(lspace, lspace, 'k--')
    ax.plot(subdata['Delta_exp'], subdata['Delta_computed'], 'o', color=color)
    
    error = subdata['Delta_exp'] - subdata['Delta_computed']
    
    ax.text(*position, '{:.2f} $\\pm$ {:.2f} ({}) '.format(numpy.mean(error), numpy.std(error), label), color=color)
        
parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-r', '--ref', default='../data/Data_XPS_C185_exp.csv')
parser.add_argument('-a', '--atom', default='C')
parser.add_argument('-o', '--output', default='Data_XPS_C185.pdf')

args = parser.parse_args()

if len(args.inputs) != len(args.names):
    raise Exception('len({}) and len({}) do not match'.format(repr(args.inputs), repr(args.names)))

data_ref = pandas.read_csv(args.ref)

data = []
graphs = []

for inp_, name in zip(args.inputs, args.names):
    data.append(prepare_data(pandas.read_csv(inp_), data_ref))
    
    g = name.split('/')[0]
    if g not in graphs:
        graphs.append(g)


figure = plt.figure(figsize=(len(graphs) * 4, 4))
axes = figure.subplots(1, len(graphs), sharey=True)

lspace = numpy.linspace(-5, 10, 101)

[ax.plot(lspace, lspace, 'k--') for ax in axes]
[ax.set_xlim(-5, 10) for ax in axes]
[ax.set_ylim(-5, 10) for ax in axes]

indexes = dict((n, 0) for n in graphs)

COLORS = ['tab:blue', 'tab:pink', 'tab:green', 'tab:red', 'tab:cyan']
POSITIONS = [(-3, x) for x in numpy.linspace(9, 5, 6)]

for data_, name in zip(data, args.names):
    graph, ref = name.split('/')
    plot_atom(axes[graphs.index(graph)], data_, 'C', COLORS[indexes[graph]], POSITIONS[indexes[graph]], ref)
    indexes[graph] += 1

[ax.set_xlabel('Experimental $\\Delta$BE (eV)') for ax in axes]
axes[0].set_ylabel('Computed $\\Delta$BE (eV)')

plt.tight_layout()
figure.savefig(args.output)

