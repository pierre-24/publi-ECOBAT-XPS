import numpy
import pathlib
import pandas

def gaussian(x, mu, sigma=0.2):
    return 1/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-.5*((x-mu)/sigma)**2)


def create_spectrum_BE(data: pandas.DataFrame, x: list, FWHM: float = 0.5, label_be: str = 'Delta_computed', label_n : str = 'N_Atoms'):
    yspace = numpy.zeros(x.shape)
    sigma = FWHM / (2*numpy.sqrt(2*numpy.log(2)))
    
    N = 0
    for BE, n in zip(data[label_be], data[label_n]):
        yspace += gaussian(x, BE, sigma) * n
        N += n
    
    if N == 0:
        return yspace
    else:
        return yspace / yspace.max()
    
def get_annotations(inp: str):
    annotations = {}
    for annotation in inp.split(','):
        if annotation == '':
            continue
        
        ax = annotation.split('=')
        df = ax[0]
        label = '='.join(ax[1:])
        atom, system = df.split('@')
        symbol, _ = atom.split('_')
        
        if symbol not in annotations:
            annotations[symbol] = {}
        
        if system not in annotations[symbol]:
            annotations[symbol][system] = []
            
        annotations[symbol][system].append((atom, label))
    
    return annotations

def annotate_graph(ax, data: pandas.DataFrame, annotations: list, x: float, y: float, color='black', position: str = 'top', shift=.1, fontsize=8):
    for a_atom, a_label in annotations:
        d = data[data['Atom_indices'].str.contains(a_atom)]
        if len(d) > 0:
            v = d.iloc[0]['Delta_computed']
            iloc = int((v - x[0]) / (x[-1] - x[0]) * len(x))
            ax.text(
                v, 
                y[iloc] + (shift if position == 'top' else -shift),
                '{lab}\n{val:.1f}' .format(lab=a_label, val=v), 
                va=('bottom' if position == 'top' else 'top'), 
                ha='center', 
                color=color,
                # rotation=90,
                fontsize=fontsize
            )


