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
        

def plot_atom(ax, data: pandas.DataFrame, atom: str):
    subdata = data[data['Atom'].str.contains(atom)]
    
    lspace = numpy.linspace(-5, 10, 101)
    ax.plot(lspace, lspace, 'k--')
    ax.plot(subdata['Delta_exp'], subdata['Delta_computed'], 'o', label='SJ')
    
    error = subdata['Delta_exp'] - subdata['Delta_computed']
    
    ax.text(-4, 8, '{} 1s (N={})'.format(atom, len(error)), fontsize=12)
    ax.text(-2, 6, '{:.2f} $\\pm$ {:.2f}'.format(numpy.mean(error), numpy.std(error)), color='tab:blue')
    
    ax.set_xlim(-5, 10)
        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_XPS_C185_SJn_0.csv')
parser.add_argument('-r', '--ref', default='../data/Data_XPS_C185_exp.csv')
parser.add_argument('-o', '--output', default='Data_XPS_C185.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)
data_ref = pandas.read_csv(args.ref)

subdata = prepare_data(data, data_ref)

figure = plt.figure(figsize=(7, 10))
(ax1, ax2), (ax3, ax4), (ax5, ax6) = figure.subplots(3, 2, sharey=True)
figure.delaxes(ax6)

ax1.set_ylim(-5, 10)

plot_atom(ax1, subdata, 'C')
plot_atom(ax2, subdata, 'N')
plot_atom(ax3, subdata, 'O')
plot_atom(ax4, subdata, 'B')
plot_atom(ax5, subdata, 'F')

[ax.set_ylabel('Computed $\\Delta$BE (eV)') for ax in [ax1, ax3, ax5]]
[ax.set_xlabel('Experimental $\\Delta$BE (eV)') for ax in [ax4, ax5]]

plt.tight_layout()
figure.savefig(args.output)

