import pandas
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

TO_J_M2 = 1.602177E-019/(0.0000000001)**2

def make_surface_energy_plot(ax, data: pandas.DataFrame, compound: str, orientations: list, areas: list):
    """Fit `E(N) = N*E + 2 * gamma * A`, where `gamma` is the surface energy.
    """
    for orientation, area in zip(orientations, areas):
        subdata = data[data['Compound'] == '{}/{}'.format(compound, orientation)]
        y = []
        for i in range(1, subdata.shape[0]):
            sdata = subdata.iloc[:i+1]
            params, _ = curve_fit(lambda x, a, b: x*a+b, sdata['N'], sdata['Energy'], ftol=1e-10)
            y.append(params[1] / (2 * area) * TO_J_M2)
        
        print('{}/{} {:.3f}'.format(compound, orientation, y[-1]))
        ax.plot(subdata.iloc[1:]['N'], y, 'o-', label='({})'.format(orientation.replace(',', '')))
    
    ax.set_title(compound)
    ax.set_xlabel("N")
    ax.legend()
        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='../data/Data_surface_energy.csv')
parser.add_argument('-o', '--output', default='Data_surface_energy.pdf')

args = parser.parse_args()

data = pandas.read_csv(args.input)

figure = plt.figure(figsize=(10,5))
ax1, ax2, ax3 = figure.subplots(1, 3, sharey=True)
make_surface_energy_plot(ax1, data, 'Ca', ['100', '110', '111'], [15.1915, 21.484, 13.1562])
ax1.set_ylabel("Surface energy (J/m$^2$)")
ax1.set_ylim(0.4, 1.25)

make_surface_energy_plot(ax2, data, 'CaO', ['100'], [23.3344, 16.4999, 10.1041])

make_surface_energy_plot(ax3, data, 'CaH2', ['100', '110', '111'], [39.9364,39.9364,39.9364])
ax3.set_title('CaH$_2$')

plt.tight_layout()
figure.savefig(args.output)
