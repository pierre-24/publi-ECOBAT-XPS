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

class Annotation:
    def __init__(self, label, x, y, sx: float = .0, sy: float = 5):
        self.label = label
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy
        self.arrow = False
    
    def annotate(self, ax, color='black', position: str = 'top', fontsize=8):
        ax.annotate(
            '{lab}\n{val:.1f}' .format(lab=self.label, val=self.x), 
            (self.x, self.y), 
            (self.sx, self.sy if (position == 'top' or self.sx != 0) else -self.sy),
            textcoords='offset pixels',
            va=('bottom' if position == 'top' else 'top'), 
            ha='center', 
            color=color,
            arrowprops=dict(facecolor='black', headwidth=3, headlength=3, width=1, color=color),
            fontsize=fontsize
        )

def annotate_graph(ax, data: pandas.DataFrame, annotations: list, x: float, y: float, color='black', position: str = 'top', shift=5, fontsize=8, mindx=.2):
    to_annotate = []
    for a_atom, a_label in annotations:
        d = data[data['Atom_indices'].str.contains(a_atom)]
        if len(d) > 0:
            v = d.iloc[0]['Delta_computed']
            iloc = int((v - x[0]) / (x[-1] - x[0]) * len(x))
            to_annotate.append(Annotation(a_label, v, y[iloc]))
    
    _annotate_graph(ax, to_annotate, color, position, shift, fontsize, mindx)

def annotate_graph2(ax, annotations: list, x: float, y: float, color='black', position: str = 'top', shift=5, fontsize=8, mindx=.2):
    """Same as annotate_graph, but using the mean of a bunch of data"""
    
    to_annotate = []
    for a_data, a_label in annotations:
        v = a_data['Delta_computed'].mean()
        iloc = int((v - x[0]) / (x[-1] - x[0]) * len(x))
        to_annotate.append(Annotation(a_label, v, y[iloc]))
    
    _annotate_graph(ax, to_annotate, color, position, shift, fontsize, mindx)
            

def _annotate_graph(ax, to_annotate: list, color='black', position: str = 'top', shift=5, fontsize=8, mindx=.2):
    to_annotate.sort(key=lambda x: x.x)
    for i in range(len(to_annotate)):
        if i == 0:
            continue
        
        if (to_annotate[i].x - to_annotate[i-1].x) < mindx:
            to_annotate[i-1].sx += 15
            to_annotate[i-1].arrow = True
            to_annotate[i].sx -= 15
            to_annotate[i].arrow = True
            
    for annotation in to_annotate:
        annotation.annotate(ax, color, position, fontsize)


