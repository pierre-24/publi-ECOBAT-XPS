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
        return yspace / N


