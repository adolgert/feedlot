'''
This creates colormaps for Matplotlib, or manipulates
the ones it provides.

A Matplotlib colors.Colormap is a callable object that
converts values normalized to (0,1) into either rgb or rgba,
which are either byte values or themselves floats in (0,1).

The colors.Colormap depends on a lookup table (_lut), that
is defined in subclasses, such as colors.LinearSegmentedColormap.
This lookup table isn't typically initialized until plot() is called.

There are two separate ways to specify a LinearSegmentedColormap.
Either give it line segments to describe red, green, and blue,
or give it functions that take (0,1) and return (0,1).
'''

import logging
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm


logger=logging.getLogger('custom_colormap')

def reduce_color_map(color_map, M):
    # Reduce the size of the colormap.
    color_map.N=M
    color_map._i_under=M
    color_map._i_over=M+1
    color_map._i_bad=M+2
    logger.debug('cmap lut shape %s' % str(color_map))
    logger.debug(color_map._lut.shape)
    return color_map


def helixmap():
    # Construct our own cubehelix color map.
    from matplotlib._cm import cubehelix
    ch_data=cubehelix(1.0, 0.5, -1.5, 0.0)
    ch_data=cm.revcmap(ch_data)
    helix_map=colors.LinearSegmentedColormap('helix', ch_data, 256)
    return helix_map



def dark_gray_cmap(m=0.3, reversed=False):
    '''
    Argument is the maximum brightness.
    '''
    _gray_data = {'red':   ((0., 0., 0.),
                            (1.0, m, m)),
              'green': ((0., 0., 0.),
                        (1.0, m, m)),
              'blue':  ((0., 0., 0.),
                        (1.0, m, m))}
    if reversed:
        _gray_data=cm.revcmap(_gray_data)
    gray_cmap=colors.LinearSegmentedColormap('lightgray', _gray_data, 256)
    return gray_cmap



def light_gray_cmap(m=0.3, gamma=1.0, reversed=False):
    '''
    This maps 0 to white (1.0) and 1 to the gray level
    specified by m. The gamma determines whether to darken
    the lower or upper part of the spectrum.
    '''
    def single(x):
        return 1 - m*(x**gamma)

    _gray_data = {'red':   single,
              'green': single,
              'blue':  single}
    if reversed:
        _gray_data=cm.revcmap(_gray_data)
    gray_cmap=colors.LinearSegmentedColormap('lightgray', _gray_data, 256)
    return gray_cmap
