import pandas
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

TO_J_M2 = 1.602177E-019/(0.0000000001)**2

def plot_ENCUT(ax, data: pandas.DataFrame, atom: str, label: str):
    ax.plot(data['ENCUT'], -data['E1s {}'.format(atom)], 'o-')
    
    ax.set_xlabel('ENCUT (eV)')
    ax.set_ylabel('Absolute BE of {}'.format(atom))
    ax.set_title(label)

def plot_kp(ax, data: pandas.DataFrame, atom: str, label: str):
    ax.plot(data['imesh'], -data['E1s {}'.format(atom)], 'o-')
    
    ax.set_xlabel('k-point mesh grid')
    ax.set_ylabel('Absolute BE of {}'.format(atom))
    
    ax.set_xticks(data['imesh'], data['mesh'])
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_title(label)
        
parser = argparse.ArgumentParser()
parser.add_argument('--input-encut', default='../data/Data_ENCUT.csv')
parser.add_argument('--input-kp', default='../data/Data_kp.csv')
parser.add_argument('-o', '--output', default='Data_conv.pdf')

args = parser.parse_args()

data_encut = pandas.read_csv(args.input_encut)
data_kp = pandas.read_csv(args.input_kp)

figure = plt.figure(figsize=(8,8))
(ax1, ax2), (ax3, ax4) = figure.subplots(2, 2)

plot_ENCUT(ax1, data_encut, 'C', 'Carbon in methanol') 
plot_ENCUT(ax2, data_encut, 'O', 'Oxygen in methanol') 
plot_kp(ax3, data_kp, 'Ca', 'Calcium in Ca(100) 1x1 slab with 6 layers') 

figure.delaxes(ax4)

plt.tight_layout()
figure.savefig(args.output)
