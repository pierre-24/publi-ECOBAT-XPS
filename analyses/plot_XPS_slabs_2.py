import pandas
import matplotlib.pyplot as plt
import numpy
import argparse

from XPS.commons import create_spectrum_BE, get_annotations, annotate_graph


def prepare_data(data: pandas.DataFrame, data_height: pandas.DataFrame, systems: list):
    REFS = {  # references
        'Ca': data[(data['System'] == 'CaO_slab/8') & (data['Atom_indices' ].str.contains('Ca_025'))].iloc[0],
        'O':  data[(data['System'] == 'ref_O')].iloc[0]
    }
    
    delta_computed = []
    
    data_out = data[data['System'].isin(systems)]
    
    for line in data_out.itertuples():
        atom = line.Atom
        delta_computed.append(line.Value - REFS[atom]['Value'])
    
    data_out.insert(5, 'Delta_computed', delta_computed)
    
    return data_out
    

def plot_atom(ax, data: pandas.DataFrame, systems: list, atom: str, color: str, label: str, shift_y: float, xrange: tuple, annotations: dict):
    
    subdata_sys1 = data[(data['System'] == systems[0]) & (data['Atom'] == atom)]
    subdata_sys2 = data[(data['System'] == systems[1]) & (data['Atom'] == atom)]
        
    lspace = numpy.linspace(*xrange, 200)
    ax.plot(xrange, [shift_y, shift_y], '-', color='grey')

    y_sys1 = create_spectrum_BE(subdata_sys1, lspace)
    ax.plot(lspace, y_sys1 + shift_y, color=color, label=label)
    y_sys2 = create_spectrum_BE(subdata_sys2, lspace)
    ax.plot(lspace, shift_y - y_sys2, '--', color=color)
    
    ax.set_xlim(*xrange)
    
    if systems[0] in annotations:
        annotate_graph(ax, subdata_sys1, annotations[systems[0]], lspace, shift_y + y_sys1, color=color)
    if systems[1] in annotations:
        annotate_graph(ax, subdata_sys2, annotations[systems[1]], lspace, shift_y - y_sys2, color=color, position='bottom')
        
parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
parser.add_argument('--height', default='../data/results_slabs/slabs_height.csv')
parser.add_argument('-s', '--systems', nargs='*', default=['CaO_OH2_slab/3', 'CaO_slab/3'])
parser.add_argument('-a', '--atoms', nargs='*')
parser.add_argument('-n', '--names', nargs='*', required=True)
parser.add_argument('-o', '--output', default='Data_XPS_slabs.pdf')
parser.add_argument('-x', '--annotate', type=get_annotations)

args = parser.parse_args()

if len(args.inputs) != len(args.names):
    raise Exception('len({}) and len({}) do not match'.format(repr(args.inputs), repr(args.names)))

data_height = pandas.read_csv(args.height)

data = []
for inp in args.inputs:
    data.append(prepare_data(pandas.read_csv(inp), data_height, args.systems))

figure = plt.figure(figsize=(len(args.atoms) * 5, len(args.names) * 1.4))
axes = figure.subplots(1, len(args.atoms), sharey=True)

COLORS = ['tab:blue', 'tab:pink', 'tab:green', 'tab:red', 'tab:cyan']

YS = 4

for i, subdata in enumerate(data):
    for j, atom_d in enumerate(args.atoms):
        a, ami, amax = atom_d.split(':')
        plot_atom(axes[j], subdata, args.systems, a, COLORS[i], args.names[i], -i * YS, (float(ami), float(amax)), annotations = args.annotate[a])
        
        if i == 0:
            axes[j].text(.05, .95, '{} {}'.format(a, '2s' if a == 'Ca' else '1s'), transform=axes[j].transAxes, fontsize=14, ha='left')

axes[1].legend()
axes[0].set_ylim(- len(args.names) * YS + 1.5, 3)
[ax.invert_xaxis() for ax in axes]
[ax.tick_params('y', left=False, labelleft=False) for ax in axes]
[ax.set_xlabel('Computed $\\Delta$BE (eV)') for ax in axes]

plt.tight_layout()
figure.savefig(args.output)

