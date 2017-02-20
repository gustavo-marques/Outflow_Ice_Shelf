# All the plotting related functions used in the notebooks are kept here
# author: Gustavo M. Marques

from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append('~/python/pyGVtools/')
import m6toolbox

# Define a function to plot a section
def plot_section(file_handle, record, xq, i=2, variable='temp',eta='e',yvar='lath',clim=(-1,0), plot_grid=True, rep='pcm', xlim=(0,200), ylim=(-1000,0), cmap=plt.cm.jet):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.
    
    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

    e = file_handle.variables[eta][record,:,:,i] # Vertical grid positions
    s = file_handle.variables[variable][record,:,:,i] # Scalar field to color
    y = file_handle.variables[yvar][:]
    x,z,q = m6toolbox.section2quadmesh(xq, e, s, representation=rep) # This yields three areas at twice the model resolution
    cs = plt.pcolormesh(x, z, q, cmap=cmap);
    #plt.colorbar()
    cs.set_clim(clim)
    if plot_grid:
       plt.plot(x, z.T, 'k', lw=0.2, hold=True)
       plt.plot(y,e[0,:],'k',lw=1); plt.plot(y,e[-1,:],'k',lw=1);
    if variable == 'v':
       plt.text(5,-500, r'max$(|v|)$:'+str(np.max(np.abs(s))), fontdict=font)
    plt.ylim(ylim)
    plt.xlim(xlim)
    return cs

