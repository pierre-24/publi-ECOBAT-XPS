import pandas
import matplotlib.pyplot as plt
import numpy
import argparse


def prepare_data(data: pandas.DataFrame, data_height: pandas.DataFrame):
    REFS = {  # references
        'Ca': data[(data['System'] == 'Ca_slab/8') & (data['Atom_indices' ].str.contains('Ca_037'))].iloc[0],
        'O':  data[(data['System'] == 'ref_O')].iloc[0]
    }
    
    delta_computed = []
    is_bulk = []
    is_surf = []
    
    data_out = data[data['System'] != 'ref_O']
    
    h_maxes = {}
    
    for line in data_out.itertuples():
        atom = line.Atom
        delta_computed.append(line.Value - REFS[atom]['Value'])
        
        if line.System not in h_maxes:
            h_maxes[line.System] = {}
        
        if line.Atom not in h_maxes[line.System]:
            h_maxes[line.System][line.Atom] = data_height[(data_height['System'] == line.System) & (data_height['Atom' ].str.contains(line.Atom))].max()['z_depth']
        
        isb, iss = False, False
        for a in line.Atom_indices.split(';'):
            h = data_height[(data_height['System'] == line.System) & (data_height['Atom' ] == a)].iloc[0]
            if h['z_depth'] <= h_maxes[line.System][line.Atom]:
                isb = True
            if h['z_depth'] >= h_maxes[line.System][line.Atom]:
                iss = True
                
            if iss and isb:
                break
        
        is_bulk.append(isb)
        is_surf.append(iss)
    
    data_out.insert(5, 'Delta_computed', delta_computed)
    data_out.insert(6, 'Is_bulk', is_bulk)
    data_out.insert(7, 'Is_surf', is_surf)
    
    return data_out
        

def plot_atom(ax, data: pandas.DataFrame, slab: str, atom: str, color: str, label: str):
    
    bulk_vals = []
    surf_vals = []

    for i in range(3, 9):
        subdata = data[(data['System'] == '{}_slab/{}'.format(slab, i)) & (data['Atom'] == atom) & (data['Delta_computed'] < 10) & (data['Delta_computed'] > -10)]
        surf_data = subdata[subdata['Is_surf'] == True]
        surf_vals.append((numpy.mean(surf_data['Delta_computed']), numpy.std(surf_data['Delta_computed'])))
        bulk_data = subdata[subdata['Is_bulk'] == True]
        bulk_vals.append((numpy.mean(bulk_data['Delta_computed']), numpy.std(bulk_data['Delta_computed'])))
    
    ax.errorbar(numpy.arange(3, 9) * (4 if slab == 'CaH2' else 2) - .1, [x[0] for x in bulk_vals], fmt='o-', color=color, yerr=[x[1] for x in bulk_vals], label=label)
    ax.errorbar(numpy.arange(3, 9) * (4 if slab == 'CaH2' else 2) + .1, [x[0] for x in surf_vals], fmt='o--', fillstyle='none', color=color, yerr=[x[1] for x in surf_vals])
        
parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('--height', default='../data/results_slabs/slabs_height.csv')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')

args = parser.parse_args()

if len(args.inputs) != len(args.names):
    raise Exception('len({}) and len({}) do not match'.format(repr(args.inputs), repr(args.names)))


data_height = pandas.read_csv(args.height)

data = []
for inp in args.inputs:
    data.append(prepare_data(pandas.read_csv(inp), data_height))

figure = plt.figure(figsize=(10, 10))
axes = figure.subplots(2, 2, sharey=True)

COLORS = ['tab:blue', 'tab:pink', 'tab:green', 'tab:red', 'tab:cyan']

for i, subdata in enumerate(data):
    plot_atom(axes[0, 0], subdata, 'Ca', 'Ca', COLORS[i], args.names[i])
    plot_atom(axes[0, 1], subdata, 'CaH2', 'Ca', COLORS[i], args.names[i])
    plot_atom(axes[1, 0], subdata, 'CaO', 'Ca', COLORS[i], args.names[i])
    plot_atom(axes[1, 1], subdata, 'CaO', 'O', COLORS[i], args.names[i])

[ax.legend() for ax in axes.flatten()]
[ax.set_ylabel('Computed $\\Delta$BE (eV)') for ax in [axes[0, 0], axes[1, 0]]]
[ax.set_xlabel('$N$') for ax in [axes[0, 1], axes[1, 1]]]

[ax.text(.05, .95, '{} / {}'.format(slab, atom), transform=ax.transAxes, fontsize=14) for ax, slab, atom in [
    (axes[0, 0], 'Ca', 'Ca 2s'),
    (axes[0, 1], 'CaH$_2$', 'Ca 2s'),
    (axes[1, 0], 'CaO', 'Ca 2s'),
    (axes[1, 1], 'CaO', 'O 1s'),
]]

plt.tight_layout()
figure.savefig(args.output)

