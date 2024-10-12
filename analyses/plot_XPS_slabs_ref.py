import pandas
import matplotlib.pyplot as plt
import numpy
import argparse
import pathlib


def prepare_data(data: pandas.DataFrame, data_height: pandas.DataFrame):
    REFS = {  # references
        'Ca': data[(data['System'] == 'CaO_slab/8') & (data['Atom_indices' ].str.contains('Ca_025'))].iloc[0],
        'O':  data[(data['System'] == 'ref_O')].iloc[0]
    }
    
    delta_computed = []
    is_bulk = []
    is_surf = []
    is_OH = []
    
    data_out = data[data['System'] != 'ref_O']
    
    for line in data_out.itertuples():
        atom = line.Atom
        delta_computed.append(line.Value - REFS[atom]['Value'])
        has_OH2 = 'OH2' in line.System
        
        isb, iss, iso = False, False, False
        for a in line.Atom_indices.split(';'):
            h = data_height[(data_height['System'] == line.System) & (data_height['Atom' ] == a)].iloc[0]
            if has_OH2:
                if h['z_depth'] <= 0.35:
                    isb = True
                if 0.35 < h['z_depth'] <= 0.4:
                    iss = True
                if h['z_depth'] > 0.4:
                    iso = True
            else:
                if h['z_depth'] <= 0.4:
                    isb = True
                if h['z_depth'] > 0.4:
                    iss = True 
                
            if iss and isb and iso:
                break
        
        is_bulk.append(isb)
        is_surf.append(iss)
        is_OH.append(iso)
    
    data_out.insert(5, 'Delta_computed', delta_computed)
    data_out.insert(6, 'Is_bulk', is_bulk)
    data_out.insert(7, 'Is_surf', is_surf)
    data_out.insert(7, 'Is_OH', is_OH)
    
    return data_out


def plot_atom(ax, data: pandas.DataFrame, slab: str, atom: str, color: str, label: str):
    
    bulk_vals = []
    surf_vals = []

    for i in range(3, 9):
        subdata = data[(data['System'] == '{}_slab/{}'.format(slab, i)) & (data['Atom'] == atom)]
        surf_data = subdata[subdata['Is_surf'] == True]
        surf_vals.append((numpy.mean(surf_data['Delta_computed']), numpy.std(surf_data['Delta_computed'])))
        bulk_data = subdata[subdata['Is_bulk'] == True]
        bulk_vals.append((numpy.mean(bulk_data['Delta_computed']), numpy.std(bulk_data['Delta_computed'])))
    
    ax.errorbar(numpy.arange(3, 9) * (4 if slab == 'CaH2' else 2) - .1, [x[0] for x in bulk_vals], fmt='o-', color=color, yerr=[x[1] for x in bulk_vals], label=label)
    ax.errorbar(numpy.arange(3, 9) * (4 if slab == 'CaH2' else 2) + .1, [x[0] for x in surf_vals], fmt='o--', fillstyle='none', color=color, yerr=[x[1] for x in surf_vals])

def make_table(f, data: pandas.DataFrame, label: str):
    def m(dt, has_OH=False):
        bulk_data = dt[dt['Is_bulk'] == True]
        f.write('{:.2f} $\\pm$ {:.2f} & '.format(numpy.mean(bulk_data['Delta_computed']), numpy.std(bulk_data['Delta_computed'])))
        surf_data = dt[dt['Is_surf'] == True]
        f.write('{:.2f} $\\pm$ {:.2f}'.format(numpy.mean(surf_data['Delta_computed']), numpy.std(surf_data['Delta_computed'])))
        if has_OH:
            hyd_data = dt[dt['Is_OH'] == True]
            f.write(' & {:.2f} $\\pm$ {:.2f}'.format(numpy.mean(hyd_data['Delta_computed']), numpy.std(hyd_data['Delta_computed'])))
    
    f.write('\\midrule\n')
    f.write('\\multicolumn{{7}}{{c}}{{{}}}  \\\\\n'.format(label))
    for slab in ['Ca', 'CaH2', 'CaO', 'CaO_OH2']:
        subdata = data[(data['System'] == '{}_slab/3'.format(slab))]
        
        if slab == 'Ca':
            f.write('\\ce{Ca^0} &')
        elif slab == 'CaO_OH2':
            f.write('\\ce{CaO.H2O} & ')
        elif slab == 'CaH2':
            f.write('\\ce{CaH2} & ')
        else:
            f.write('{} & '.format(slab))
        
        data_Ca = subdata[subdata['Atom'] == 'Ca']
        m(data_Ca)
        
        f.write(' && ')
        
        if 'O' in slab:
            data_O = subdata[subdata['Atom'] == 'O']
            m(data_O, has_OH='OH2' in slab)
            if 'OH2' not in slab:
                f.write(' & ---')
        else:
            f.write('--- & --- & ---')
        
        f.write('\\\\\n')


parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('--height', default='../data/results_slabs/slabs_height.csv')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')
parser.add_argument('-t', '--table')

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

if args.table:
    with pathlib.Path(args.table).open('w') as f:
        for i, subdata in enumerate(data):
            make_table(f, subdata,  args.names[i])
