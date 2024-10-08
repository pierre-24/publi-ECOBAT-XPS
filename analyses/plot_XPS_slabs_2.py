import pandas
import matplotlib.pyplot as plt
import numpy
import argparse

from XPS.commons import create_spectrum_BE


def prepare_data(data: pandas.DataFrame, data_height: pandas.DataFrame):
    REFS = {  # references
        'Ca': data[(data['System'] == 'Ca_slab/8') & (data['Atom_indices' ].str.contains('Ca_037'))].iloc[0],
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
        
        isb, iss, iso = False, False, False
        for a in line.Atom_indices.split(';'):
            h = data_height[(data_height['System'] == line.System) & (data_height['Atom' ] == a)].iloc[0]
            if h['z_depth'] <= 0.35:
                isb = True
            if 0.35 < h['z_depth'] <= 0.4:
                iss = True
            if h['z_depth'] > 0.4:
                iso = True
                
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
        

def plot_atom(ax, data: pandas.DataFrame, system: str, atom: str, color: str, label: str, shift_y: float = .0, space: tuple = (-10, 10)):
    subdata = data[(data['System'] == system) & (data['Atom'] == atom)]
    
    lspace = numpy.linspace(*space, 200)
    y = create_spectrum_BE(subdata, lspace)
    ax.plot(lspace, y + shift_y, label=label, color=color)
    
    m_bulk = subdata[subdata['Is_bulk'] == True]['Delta_computed'].mean()
    ax.plot([m_bulk, m_bulk], [shift_y, shift_y + 0.3], '-', color=color)
    ax.text(m_bulk, shift_y + .3, 'b', va='bottom', ha='center', color=color)
    m_surf = subdata[subdata['Is_surf'] == True]['Delta_computed'].mean()
    ax.plot([m_surf, m_surf], [shift_y, shift_y + 0.3], '-', color=color)
    ax.text(m_surf, shift_y + .3, 's', va='bottom', ha='center', color=color)
    
    if atom == 'O':
        m_OH = subdata[subdata['Is_OH'] == True]['Delta_computed'].mean()
        ax.plot([m_OH, m_OH], [shift_y, shift_y + 0.3], '-', color=color)
        
parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('--height', default='../data/results_slabs/slabs_height.csv')
parser.add_argument('-s', '--system', default='CaO_OH2_slab/3')
parser.add_argument('-a', '--atoms', nargs='*')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')

args = parser.parse_args()

if len(args.inputs) != len(args.names):
    raise Exception('len({}) and len({}) do not match'.format(repr(args.inputs), repr(args.names)))


data_height = pandas.read_csv(args.height)

data = []
for inp in args.inputs:
    data.append(prepare_data(pandas.read_csv(inp), data_height))

figure = plt.figure(figsize=(len(args.atoms) * 5, 6))
axes = figure.subplots(1, len(args.atoms))

COLORS = ['tab:blue', 'tab:pink', 'tab:green', 'tab:red', 'tab:cyan']

YS = 1.5

for i, subdata in enumerate(data):
    for j, atom in enumerate(args.atoms):
        a, ami, amax = atom.split(':')
        plot_atom(axes[j], subdata, args.system, a, COLORS[i], args.names[i], shift_y = i * YS, space = (float(ami), float(amax)))
        
        if i == 0:
            axes[j].text(.95, .95, 'CaO$\\cdot$H$_2$O / {} {}'.format(a, '2s' if a == 'Ca' else '1s'), transform=axes[j].transAxes, fontsize=14, ha='right')

[ax.legend(loc='upper left') for ax in axes]
[ax.invert_xaxis() for ax in axes]
[ax.tick_params('y', left=False, labelleft=False) for ax in axes]
[ax.set_xlabel('Computed $\\Delta$BE (eV)') for ax in axes]

plt.tight_layout()
figure.savefig(args.output)

