import pandas
import matplotlib.pyplot as plt
import numpy
import argparse


def prepare_data(data: pandas.DataFrame, data_height: pandas.DataFrame):
    REFS = {  # references
        'Ca': data[(data['System'] == 'CaO_slab/8') & (data['Atom_indices' ].str.contains('Ca_025'))].iloc[0],
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
            if abs(h['z_depth'] - h_maxes[line.System][line.Atom]) < .05:
                iss = True
            else:
                isb = True
                
            if iss and isb:
                break
        
        is_bulk.append(isb)
        is_surf.append(iss)
    
    data_out.insert(5, 'Delta_computed', delta_computed)
    data_out.insert(6, 'Is_bulk', is_bulk)
    data_out.insert(7, 'Is_surf', is_surf)
    
    return data_out
        

def plot_atom(ax, data: pandas.DataFrame, slab: str, atom: str, color: str, marker: str, label: str, maxstd: float = .5):
    
    bulk_vals = []
    surf_vals = []

    for i in range(3, 9):
        subdata = data[(data['System'] == '{}_slab/{}'.format(slab, i)) & (data['Atom'] == atom)]
        surf_data = subdata[subdata['Is_surf'] == True]
        if numpy.std(surf_data['Delta_computed']) < maxstd:
            surf_vals.append((i, numpy.mean(surf_data['Delta_computed']), numpy.std(surf_data['Delta_computed'])))
        bulk_data = subdata[subdata['Is_bulk'] == True]
        if numpy.std(bulk_data['Delta_computed']) < maxstd:
            bulk_vals.append((i, numpy.mean(bulk_data['Delta_computed']), numpy.std(bulk_data['Delta_computed'])))
    
    ax.errorbar([x[0] - .1 for x in bulk_vals], [x[1] for x in bulk_vals], fmt=marker + '-', color=color, yerr=[x[2] for x in bulk_vals], label=label)
    ax.errorbar([x[0] + .1 for x in surf_vals], [x[1] for x in surf_vals], fmt=marker + '--', fillstyle='none', color=color, yerr=[x[2] for x in surf_vals])
        
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

NY = int(numpy.ceil(len(args.inputs) / 2))
figure = plt.figure(figsize=(8, NY * 3))
axes = figure.subplots(NY, 2, sharey=True)

COLORS = ['tab:blue', 'tab:pink', 'tab:green', 'tab:red', 'tab:cyan']

for i, (ax, subdata) in enumerate(zip(axes.flatten(), data)):
    plot_atom(ax, subdata, 'Ca', 'Ca', COLORS[i], 'o', 'Ca$^0$')
    plot_atom(ax, subdata, 'CaH2', 'Ca', COLORS[i], 's', 'CaH$_2$')
    plot_atom(ax, subdata, 'CaO', 'Ca', COLORS[i], '^', 'CaO')
    
    ax.text(.05, .9, args.names[i], transform=ax.transAxes, fontsize=12, color=COLORS[i])

[ax.legend() for ax in axes.flatten()]
[ax.set_ylabel('Computed $\\Delta$BE (eV)') for ax in [axes[0, 0], axes[1, 0]]]
[ax.set_xlabel('$N / N_0$') for ax in axes.flatten()]

if len(args.inputs) % 2 != 0:
    figure.delaxes(axes.flatten()[-1])

ymi, yma = axes[0,0].get_ylim()    
axes[0,0].set_ylim(ymi, yma + .75)

plt.tight_layout()
figure.savefig(args.output)

