import numpy
import pathlib
import pandas

def gaussian(x, mu, sigma=0.2):
    return 1/(sigma*numpy.sqrt(2*numpy.pi))*numpy.exp(-.5*((x-mu)/sigma)**2)


def create_spectrum(data: pandas.DataFrame, x: list, ref: float, sigma: float = 0.2):
    yspace = numpy.zeros(x.shape)
    
    N = 0
    for BE, n in zip(data['BE'], data['N']):
        yspace += gaussian(x, BE - ref, sigma) * n
        N += n
    
    if N == 0:
        return yspace
    else:
        return yspace / N

