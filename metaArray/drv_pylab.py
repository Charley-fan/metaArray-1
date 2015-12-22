#  drv_pylab.py
#  
#  Copyright 2015 Unknown <charley@utc2d>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

# This file contain a number of drivers classes to matplotlib.

from matplotlib.pyplot import figure, tight_layout
from matplotlib import cm
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
from matplotlib import rcParams
import matplotlib.pyplot as plt

from numpy import linspace, logspace, arange
from numpy import sign
from numpy import average
from numpy import dtype
from numpy import real, imag
from numpy import abs
from numpy import angle, pi
from numpy import log10
from numpy import unwrap
from numpy import arctan2
from numpy import floor, ceil

from scipy.signal import freqz

from misc import unitPrefix, prettyunit


def plotcomplexpolar(metaAry, axis = -1, size = (10, 7.5), dpi = 75, \
                    grid = True, legend = 0, fontsize = 14, linewidth = 2.0, \
                    log_mag=False, mag_label=None, mag_unit=None, \
                    pha_label=None, pha_unit=None, degree=False):
    """
    metaArray function to do a simple 1D plot of complex array (metaAry[axis]) as magnitude and phase angle.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    """
    
    assert type(axis) is int, "Axis is not an integer: %r" % axis
    
    fontsize = float(fontsize)
    
    if legend is None:
        legend = 0
    
    if mag_label is None:
        mag_label = "Magnitude"
    
    if pha_label is None:
        pha_label = "Phase"
    
    mag = abs(metaAry.data)
    pha = angle(metaAry.data)
    
    # Load the plotting ranges and units
    x0 = metaAry.get_range(axis, 'begin')
    x1 = metaAry.get_range(axis, 'end')
    my0 = min(mag)
    my1 = max(mag)
    py0 = floor(4 * min(pha) / pi)      # Round to the nearest pi/4
    py1 = ceil(4 * max(pha) / pi)
    
    pticks = arange(py0, py1+1)          # Ticks in pi/4 interval
    
    xunit = metaAry.get_range(axis, 'unit')
    
    if mag_unit == None:
        myunit = metaAry['unit']
    else:
        myunit = str(mag_unit)
    
    if pha_unit == None:
        if degree == True:
            pyunit = 'Deg.'
        else:
            pyunit = 'Rad.'
    else:
        pyunit = pha_unit
    
    # Leave 10% margin in the y axis
    if log_mag == True:
        my0 = log10(my0)
        my1 = log10(my1)
    
    mmean = average((my0, my1))
    mreach = abs(my0-my1) / 2 / 0.9
    my0 = sign(my0-mmean) * mreach + mmean
    my1 = sign(my1-mmean) * mreach + mmean
    
    if log_mag == True:
        my0 = 10**my0
        my1 = 10**my1
    
    pmean = average((py0, py1))
    preach = abs(py0-py1) / 2 / 0.9
    py0 = sign(py0-pmean) * preach + pmean
    py1 = sign(py1-pmean) * preach + pmean
    
    my0, my1 = scale_check(my0, my1)
    py0, py1 = scale_check(py0, py1)
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    myunit, my0, my1, myscale = prettyunit(myunit, my0, my1)
    pyunit, py0, py1, pyscale = prettyunit(pyunit, py0, py1)
    
    x = metaAry.get_axis(axis=axis)
    
    if xscale != 1: 
        x *= xscale
    
    if myscale != 1: 
        mag *= myscale
    
    if pyscale != 1: 
        pha *= pyscale
    
    xlabl = lbl_repr(metaAry.get_range(axis, 'label'), xunit)
    mylabl = lbl_repr(metaAry['label'], myunit, mag_label)
    pylabl = lbl_repr(metaAry['label'], pyunit, pha_label)
    
    title = metaAry['name']
    
    # Done the preparation, do the actual plotting
    
    fig, ax1 = plt.subplots(figsize=size, dpi=dpi)
    ax2 = ax1.twinx()
    ax1.grid(True, which="both", ls="--", color='g')
    
    ######
    
    ax1.plot(x, mag, 'b-', linewidth=linewidth, label=mag_label)
        
    if degree == True:
        ax2.plot(x, pha * 180 / pi, 'r--', linewidth=linewidth, label=pha_label)
    else:
        ax2.plot(x, pha, 'r--', linewidth=linewidth, label=pha_label)
    
    ######
    
    ax1.set_ylabel(mylabl, color='b', fontsize=fontsize)
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
        tl.set_fontsize(fontsize)
    
    if log_mag == True:
        ax1.set_yscale('log', nonposy='clip')
    
    ax1.set_ylim([my0, my1])
    
    ######
    
    ax2.set_ylabel(pylabl, color='r', fontsize=fontsize)
    
    if degree == True:
        py0 *= 45
        py1 *= 45
        
        ax2.set_ylim([py0, py1])
        ax2.set_yticks(pticks*45)
    else:
        py0 = py0 * pi / 4
        py1 = py1 * pi / 4
        
        ax2.set_ylim([py0, py1])
        ax2.set_yticks(pticks * pi / 4)
        
        pticks_lbl = []
        for pt in pticks:
            if abs(pt) == 1:
                val = r'\pi/4'
                
            elif abs(pt) == 2:
                val = r'\pi/2'
                
            elif abs(pt) == 3:
                val = r'3\pi/4'
                
            elif abs(pt) == 4:
                val = r'\pi'
            
            if abs(pt) == 0:
                val = '0'
                
            elif sign(pt) == -1:
                val = '-' + val
            
            pticks_lbl.append('$' + val + '$')
        
        ax2.set_yticklabels(pticks_lbl)
    
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
        tl.set_fontsize(fontsize)
    
    ######
    
    ax1.set_title(title, fontsize=fontsize*1.3)
    
    if metaAry.get_range(axis, 'log') == True:
        ax1.set_xscale("log", nonposx='clip')
        ax2.set_xscale("log", nonposx='clip')
    
    for tl in ax1.get_xticklabels():
        tl.set_fontsize(fontsize)
    
    ax1.set_xlabel(xlabl, fontsize=fontsize)
    ax1.set_xlim([x0, x1])

    
    if legend >= 0:
        lns1, lbl1 = ax1.get_legend_handles_labels()
        lns2, lbl2 = ax2.get_legend_handles_labels()
        ax1.legend(lns1 + lns2, lbl1 + lbl2, loc=legend)
    
    ######
    
    ## fig = figure(figsize=size, dpi = dpi)
    ## host = SubplotHost(fig, 111)
    
    #fig.add_subplot(host)
    #par = host.twinx()
    
    #host.plot(x, mag, 'b-', label=mag_label)
    #par.plot(x, pha, 'r--', label=pha_label)
    
    #host.grid(grid, which="both", color='g')
    
    #host.set_xlabel(xlabl)
    #host.set_ylabel(mylabl, color='b')
    #par.set_ylabel(pylabl, color='r')
    
    #if log_mag == True:
        #host.set_yscale('log', nonposy='clip')
    
    #host.set_xlim([x0, x1])
    #host.set_ylim([my0, my1])
    #par.set_ylim([py0, py1])
    
    
    #for tl in host.get_xticklabels():
        #tl.set_fontsize(fontsize)
    
    #for tl in host.get_yticklabels():
        #tl.set_color('b')
        #tl.set_fontsize(fontsize)
    
    #for tl in par.get_yticklabels():
        #tl.set_color('r')
        #tl.set_fontsize(fontsize)
    
    ## host.set_title(title, fontsize=int(fontsize*1.3))
    #host.set_title(title, fontsize=fontsize)
    
    #if legend >= 0:
        #host.legend(loc=legend)
    
    return fig, ax1, ax2
    


def plotcomplex(metaAry, size = (10, 7.5), dpi = 75, grid = True, \
                legend = 0, fontsize = 15, real_label=None, imag_label=None):
    """
    metaArray function to do a simple 1D plot of complex array as real and imaginary parts.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    """
    
    if legend is None:
        legend = 0
    
    if real_label is None:
        real_label = "Real"
    
    if imag_label is None:
        imag_label = "Imaginary"
    
    axis = metaAry['range']
    rdata = metaAry.data.real
    idata = metaAry.data.imag
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    ry0 = min(rdata)
    ry1 = max(rdata)
    iy0 = min(idata)
    iy1 = max(idata)
    xunit = axis['unit'][0]
    ryunit = metaAry['unit']
    iyunit = metaAry['unit']
    
    # Leave 10% margin in the y axis
    rmean = average((ry0, ry1))
    rreach = abs(ry0-ry1) / 2 / 0.9
    ry0 = sign(ry0-rmean) * rreach + rmean
    ry1 = sign(ry1-rmean) * rreach + rmean
    
    imean = average((iy0, iy1))
    ireach = abs(iy0-iy1) / 2 / 0.9
    iy0 = sign(iy0-imean) * ireach + imean
    iy1 = sign(iy1-imean) * ireach + imean
    
    ry0, ry1 = scale_check(ry0, ry1)
    iy0, iy1 = scale_check(iy0, iy1)
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    ryunit, ry0, ry1, ryscale = prettyunit(ryunit, ry0, ry1)
    iyunit, iy0, iy1, iyscale = prettyunit(iyunit, iy0, iy1)
    
    if ryscale != 1:
        rdata = rdata.copy() * ryscale
        
    if iyscale != 1:
        idata = idata.copy() * iyscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    rylabl = lbl_repr(metaAry['label'], ryunit, real_label + ' part')
    iylabl = lbl_repr(metaAry['label'], iyunit, imag_label + ' part')
    
    title = metaAry['name']
    
    fig = figure(figsize=size, dpi = dpi)
    host = SubplotHost(fig, 111)
    
    fig.add_subplot(host)
    par = host.twinx()
    
    #if axis['log'][0] == False:
    #    x = linspace(x0, x1, len(metaAry))
    #else:
    #    raise NotImplemented, "Log axis is not yet implemented."
    
    x = metaAry.get_axis()
    
    host.plot(x, rdata, 'b-', label=lbl_repr(axis['label'][0], '', real_label))
    par.plot(x, idata, 'r--', label=lbl_repr(axis['label'][0], '', real_label))
    
    host.grid(grid)
    
    host.set_xlabel(xlabl, fontsize=fontsize)
    host.set_ylabel(rylabl, fontsize=fontsize)
    par.set_ylabel(iylabl, fontsize=fontsize)
    
    host.set_xlim([x0, x1])
    host.set_ylim([ry0, ry1])
    par.set_ylim([iy0, iy1])
    
    if fontsize is not None:
        host.set_title(title, fontsize=int(fontsize*1.3))
    else:
        host.set_title(title)
        
    if legend >= 0:
        host.legend(loc=legend)
    
    return fig, host, par
    



def multiplot(metaAry, size = (10, 7.5), dpi = 75, grid = True, \
                legend = 0, fontsize = 15, real_label=None, imag_label=None, \
                fig=None, host=None, par=None):
    """
    metaArray function to do a simple 1D plot of complex array as real and imaginary parts.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    """
    
    if legend is None:
        legend = 0
    
    if real_label is None:
        real_label = "Real"
    
    if imag_label is None:
        imag_label = "Imaginary"
    
    axis = metaAry['range']
    rdata = metaAry.data.real
    idata = metaAry.data.imag
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    ry0 = min(rdata)
    ry1 = max(rdata)
    iy0 = min(idata)
    iy1 = max(idata)
    xunit = axis['unit'][0]
    ryunit = metaAry['unit']
    iyunit = metaAry['unit']
    
    # Leave 10% margin in the y axis
    rmean = average((ry0, ry1))
    rreach = abs(ry0-ry1) / 2 / 0.9
    ry0 = sign(ry0-rmean) * rreach + rmean
    ry1 = sign(ry1-rmean) * rreach + rmean
    
    imean = average((iy0, iy1))
    ireach = abs(iy0-iy1) / 2 / 0.9
    iy0 = sign(iy0-imean) * ireach + imean
    iy1 = sign(iy1-imean) * ireach + imean
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    ryunit, ry0, ry1, ryscale = prettyunit(ryunit, ry0, ry1)
    iyunit, iy0, iy1, iyscale = prettyunit(iyunit, iy0, iy1)
    
    if ryscale != 1:
        rdata = rdata.copy() * ryscale
        
    if iyscale != 1:
        idata = idata.copy() * iyscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    rylabl = lbl_repr(metaAry['label'], ryunit, real_label + ' part')
    iylabl = lbl_repr(metaAry['label'], iyunit, imag_label + ' part')
    
    title = metaAry['name']
    
    fig = figure(figsize=size, dpi = dpi)
    host = SubplotHost(fig, 111)
    
    fig.add_subplot(host)
    par = host.twinx()
    
    #if axis['log'][0] == False:
    #    x = linspace(x0, x1, len(metaAry))
    #else:
    #    raise NotImplemented, "Log axis is not yet implemented."
    
    x = metaAry.get_axis()
    
    host.plot(x, rdata, 'b-', label=lbl_repr(axis['label'][0], '', real_label))
    par.plot(x, idata, 'r--', label=lbl_repr(axis['label'][0], '', real_label))
    
    host.grid(grid)
    
    host.set_xlabel(xlabl, fontsize=fontsize)
    host.set_ylabel(rylabl, fontsize=fontsize)
    par.set_ylabel(iylabl, fontsize=fontsize)
    
    host.set_xlim([x0, x1])
    host.set_ylim([ry0, ry1])
    par.set_ylim([iy0, iy1])
    
    if fontsize is not None:
        host.set_title(title, fontsize=int(fontsize*1.3))
    else:
        host.set_title(title)
        
    if legend >= 0:
        host.legend(loc=legend)
    
    return fig, host, par
    




def plot1d(metaAry, size = (10, 7.5), dpi = 75, grid = True, legend = None, fontsize = 15,\
            fig = None, ax = None, label = None):
    """
    metaArray function to do a simple 1D plot.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    
    label   Label for the legend display, default to metaAry['range']['label'][0]
    
    """
    
    if metaAry.dtype is dtype('complex'):
        return plotcomplex(metaAry, size = size, dpi = dpi, grid = grid, legend = legend, fontsize = fontsize)
        
    if legend is None:
        legend = -1
        
    axis = metaAry['range']
    data = metaAry.data
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    y0 = min(metaAry.data)
    y1 = max(metaAry.data)
    xunit = axis['unit'][0]
    yunit = metaAry['unit']
    
    # Leave 10% margin in the y axis
    mean = average((y0, y1))
    reach = abs(y0-y1) / 2 / 0.9
    y0 = sign(y0-mean) * reach + mean
    y1 = sign(y1-mean) * reach + mean
    
    y0, y1 = scale_check(y0, y1)
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    yunit, y0, y1, yscale = prettyunit(yunit, y0, y1)
    
    if yscale != 1:
        data = data.copy() * yscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    ylabl = lbl_repr(metaAry['label'], yunit)
    
    title = metaAry['name']
    
    # check if object is 1D metaArray object
    if fig is None:
        fig = figure(figsize=size, dpi = dpi)
    
    if ax is None:
        ax = fig.add_subplot(111)
    else:
        x00, x01 = ax.get_xlim()
        y00, y01 = ax.get_ylim()
        
        x0 = min((x0, x00))
        y0 = min((y0, y00))
        x1 = max((x1, x01))
        y1 = max((y1, y01))
    
    if axis['log'][0] == False:
        x = linspace(x0, x1, len(metaAry))
    else:
        raise NotImplemented
    
    if label is None:
        label = axis['label'][0]
    
    ax.plot(x, data, label=label)
    
    ax.grid(grid)
    
    ax.set_xlabel(xlabl, fontsize=fontsize)
    ax.set_ylabel(ylabl, fontsize=fontsize)
    
    ax.set_xlim([x0, x1])
    ax.set_ylim([y0, y1])
    
    if fontsize is not None:
        ax.set_title(title, fontsize=int(fontsize*1.3))
    else:
        ax.set_title(title)
        
    if legend >= 0:
        ax.legend(loc=legend)
    
    return fig, ax
    


def plot2d(metaAry, size = (10, 7.5), dpi = 75, fontsize = 15, cmap = None, \
            nticks = 5, aspect_ratio = 1.0, corient = 'vertical', cformat = None, 
            vmin = None, vmax = None, interpolation = 'sinc'):
    """
    metaArray function to do a simple 2D plot.
    
    size            Plot size, default to (10, 7.5)
    dpi             Dot Per Inch for raster graphics
    fontsize        Norminal font size
    cmap            Colour map, default is pyplot.cm.spectral
    nticks          Number of ticks in the colour bar
    aspect_ratio    Aspect ratio of the plot {float|'ij'|'xy'}
                        float:  Fixed aspect ratio by the given number
                        'ij':   Same aspect ratio as ij space
                        'xy':   Same aspect ratio as xy space 
    corient         Colorbar orientation ('vertical'|'horizontal')
    cformat         Colorbar format [ None | format string | Formatter object ]
    vmin            Minimum value for the colour scale
    vmax            Maximum value for the coloir scale
    interpolation   Colour interpolation methods 
                    [None, 'none', 'nearest', 'bilinear', 'bicubic',
                    'spline16', 'spline36', 'hanning', 'hamming', 'hermite',
                    'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 
                    'mitchell', 'sinc', 'lanczos']
    """
    
    if cmap is None:
        cmap = cm.spectral
    
    if corient is not 'horizontal':
        corient = 'vertical'
    
    axis = metaAry['range']
    data = metaAry.data
    
    x0 = axis['begin'][0] 
    x1 = axis['end'][0]
    y0 = axis['begin'][1]
    y1 = axis['end'][1]
    
    # Try to work out the aspect ratio and plot size before scailing the axis
    # Aspect ratio of the plot
    if aspect_ratio == 'ij':
        ratio = data.shape
        ratio = float(ratio[1]) / ratio[0]
    elif aspect_ratio == 'xy':
        ratio = float(y1 - y0) / (float(x1 - x0))
    else:
        try:
            ratio = float(aspect_ratio)
        except:
            print "*** Warning! Unrecognisable aspect ratio spec. Using the default instead."
            ratio = 1.0
    
    ## Make plot with vertical (default) colorbar
    #if size == 'default':
        ## Default size aspect ratio is 4:3, and size to be (10, 7.5)
        ## Adjust size according to image aspect ratio
                
        #if corient == 'vertical':
            #size_ratio = ratio * 0.75
            #size = (10, 10 * size_ratio)
        #else:
            #size_ratio = ratio / 0.75
            #size = (7.5 / size_ratio, 7.5)
        
        ##diagonal =  12.5**2
        ##size_x = (diagonal / (1 + size_ratio**2))**0.5
        ##size_y = (diagonal / size_x**2)**0.5
        
        ##size = (size_x, size_y)
    
    
    # Try to work out the colour scale
    if vmin is None:
        v0 = metaAry.data.min()
    else:
        v0 = vmin
    
    if vmax is None:
        v1 = metaAry.data.max()
    else:
        v1 = vmax
    
    # In case the values are all the same
    v0, v1 = scale_check(v0, v1)
    
    xunit = axis['unit'][0]
    yunit = axis['unit'][1]
    vunit = metaAry['unit']
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    yunit, y0, y1, yscale = prettyunit(yunit, y0, y1)
    vunit, v0, v1, vscale = prettyunit(vunit, v0, v1)
    
    if vscale != 1:
        data = data.copy() * vscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    ylabl = lbl_repr(axis['label'][1], yunit)
    vlabl = lbl_repr(metaAry['label'], vunit)
    
    ticks = linspace(v0, v1, nticks)
    ticks_lbl = []
    
        
    # Matplotlib inshow data in transposed from metaArray convension
    # And it adjust the aspect ratio based on the prefix corrected number
    ratio /= float(y1 - y0) / float(x1 - x0)
    ratio = abs(ratio)          #  This is the number fed to matplotlib
    
    for i in range(nticks):
        ticks_lbl.append("%(val)0.4g" % {'val':ticks[i]})
    
    fig = figure(figsize=size, dpi = dpi)
    ax = fig.add_subplot(111)
    
    extent = (x0, x1, y0, y1)
    cax = ax.imshow(data.transpose()[::-1], cmap=cmap, extent=extent, interpolation = interpolation, vmin = v0, vmax = v1, aspect=ratio)
    cbar = fig.colorbar(cax, ticks=ticks, orientation=corient, format=cformat)
    
    # ax.set_size(fontsize)
    ax.set_xlabel(xlabl, fontsize=fontsize)     #   Label font size
    ax.set_ylabel(ylabl, fontsize=fontsize)
    rcParams.update({'font.size': fontsize})    #   Value font size
    
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar.ax.set_yticklabels(ticks_lbl)
    cbar.set_label(vlabl, fontsize=fontsize)
    
    if fontsize is not None:
        ax.set_title(metaAry['name'], fontsize=int(fontsize*1.3))
    else:
        ax.set_title(metaAry['name'])
    
    # fig.tight_layout()
    
    return fig, ax
    

def lbl_repr(label = None, unit = None, string = None):
    """
    Format axis label and unit into a nice looking string
    
    String: Additional string between label and unit.
    
    "label [string] (unit)"
    """
    lbl = ''
    try:
        # Ignore label if it is not a string, for it can be None also
        lbl += label
    except TypeError:
        pass
        
    try:
        # Append the additional arguement if exist
        lbl += ' [' + string + ']'
    except TypeError:
        pass
    
    try:
        if unit == '':
            pass                    # Unit less quantities
        else:
            lbl += ' ('  + unit + ')'
    except TypeError:
        # Most likely unit is not defined, i.e. not a string.
        lbl += ' (Arb.)'
        
    return lbl
    


def scale_check(v0, v1):
    """
    Check if the scale limits are identical, if so, return a 0.1% difference
    between the limits
    """
    v0 = float(v0)
    v1 = float(v1)
    
    # In case the values are all the same
    if v0 == v1:
        if v0 == 0:
            v0 -= 0.0005
            v1 += 0.0005
        else:
            v0 -= v0*0.0005
            v1 += v1*0.0005
    
    return v0, v1




def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    
    return

def multiplot(metaArylst, size = (10, 7.5), dpi = 75, grid = True, \
                legend = 0, fontsize = 15, mode='inc', \
                unitChk = True, title = False, group_y = True, \
                fig = None, ax = None):
    """
    metaArray function to put the specified list of metaArrays in to a 1D plot.
    
    metaArylst      [metAry1, metAry2, ... metAryN] List of 1-D metaAry
    size            Plot size, default to (10, 7.5)
    dpi             Dot Per Inch for raster graphics
    grid            Display grid
    fontsize        Norminal font size
    mode            ['inc'|'exc'] 
    unitChk         Whether to check if all metArys has the same x-axis unit
    title           Whether to generate figure title, if title is string, use it as is.
    group_y         Identical metArys['unit'] are default to use the same y-axis scale
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    
    mode:
        'inc'       Inclusive, plot all data (x-axis) ranges
        'exc'       Exclusive, plot only common (x-axis) ranges to all data
    
    
    label       Label for the legend display, default to metaAry['range']['label'][0]
    
    """
        
    if legend is None: legend = -1
    
    # First pass, work out the overall axis range
    
    xunit = metaArylst[0].get_range(0, 'unit')
    
    x0 = metaArylst[0].get_range(0, 'begin')
    x1 = metaArylst[0].get_range(0, 'end')
    
    # unitChk Make sure they have identical x-units:
    if unitChk:
        for metAry in metaArylst:
            if metAry.get_range(0, 'unit') != xunit: 
                raise UnitError, "Axis unit description do no match (" + str(metAry['name']) + ")"
    
    # Identify the x ranges
    if mode is 'inc':
        for metAry in metaArylst:
            if metAry.get_range(0, 'begin') < x0: x0 = metAry.get_range(0, 'begin')
            if metAry.get_range(0, 'end') > x1: x1 = metAry.get_range(0, 'end')
    else:
        for metAry in metaArylst:
            if metAry.get_range(0, 'begin') > x0: x0 = metAry.get_range(0, 'begin')
            if metAry.get_range(0, 'end') < x1: x1 = metAry.get_range(0, 'end')
    
    # Define the x-axis
    # x0
    # x1
    # xunit
    # xscale
    # xlabl
    xunit, x0_scl, x1_scl, xscale = prettyunit(xunit, x0, x1)
    xlabl = lbl_repr(axis['label'][0], xunit)
    
    # Check how many y-axies are necessary
    y_lst = {}
    for idx in range(len(metaArylst)):
        
        metAry = metaArylst[idx]
        yunit = metAry['unit']
        
        if group_y:
            
            if y_lst.has_key(yunit):
                
                if y_lst[yunit]['y0'] > metAry.data.min():
                    y_lst[yunit]['y0'] = metAry.data.min()
                
                if y_lst[yunit]['y1'] < metAry.data.max():
                    y_lst[yunit]['y1'] = metAry.data.max()
            else:
                y_lst[yunit] = {}
                y_lst[yunit]['y0'] = metAry.data.min()
                y_lst[yunit]['y1'] = metAry.data.max()
                y_lst[yunit]['yunit'] = yunit
                y_lst[yunit]['ylabl'] = metaAry['label']
        else:
            idx_str = str(idx)
            y_lst[idx_str] = {}
            y_lst[idx_str]['y0'] = metAry.data.min()
            y_lst[idx_str]['y1'] = metAry.data.max()
            y_lst[idx_str]['yunit'] = yunit
            y_lst[idx_str]['ylabl'] = metaAry['label']
    
    # Per y-axis defnitions
    # y0
    # y1
    # yunit
    # yscale
    # ylabl
    for ykey in y_lst.keys():
        
        y0 = y_lst[ykey]['y0']
        y1 = y_lst[ykey]['y1']
        yunit = y_lst[ykey]['yunit']
        ylabl = y_lst[ykey]['ylabl']
        
        # Leave 10% margin in the y axis
        ymean = average((y0, y1))
        yreach = abs(y0-y1) / 2 / 0.9
        y0 = sign(y0-mean) * reach + mean
        y1 = sign(y1-mean) * reach + mean
        
        y0, y1 = scale_check(y0, y1)
        
        # Apply unit prefix if unit is defined
        yunit, y0, y1, yscale = prettyunit(yunit, y0, y1)
        
        # Update the label
        ylabl = lbl_repr(ylabl, yunit)
        
        y_lst[ykey]['y0'] = y0
        y_lst[ykey]['y1'] = y1
        y_lst[ykey]['yunit'] = yunit
        y_lst[ykey]['yscale'] = yscale
        y_lst[ykey]['ylabl'] = ylabl
    
    
    #Generate the fig
    if fig is None:
        fig = figure(figsize=size, dpi = dpi)
        fig.subplots_adjust(right=0.75)
    
    if ax is None:
        ax0 = fig.add_subplot(111)
    else:
        ax0 = ax
    
    # Generate the ax
    ykeys = y_lst.keys()
    ykey = ykeys[0]
    y_lst[ykey]['ax'] = ax0
    
    for ykey in ykeys[1:]:
        y_lst[ykey]['ax'] = ax0.twinx()
    
    for ykey in ykeys[2:]:
        ax = y_lst[ykey]['ax']
         
        # Offset the right spine of par2.  The ticks and label have already been
        # placed on the right by twinx above.
        ax.spines["right"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(ax)
        # Second, show the right spine.
        ax.spines["right"].set_visible(True)
    
    # Plot each metaArray in the list
    line_lst = []
    for idx in range(len(metaArylst)):
        # x0, x1
        # xunit, xscale, xlabl, x0_scl, x1_scl
        
        if mode is 'inc':
            metAry = metaArylst[idx]
        else: # Exclusive mode
            metAry = metaArylst[idx][float(x0):float(x1)]
        
        if group_y:
            ykey = metAry['unit']
        else:
            ykey = str(idx)
        
        # y0, y1, yscale, yunit, ylabl, ax
        yscale = y_lst[ykey]['yscale']
        ax = y_lst[ykey]['ax']
        
        # Obtain the scaled data pair
        x_axis = metaAry.get_axis() * xscale
        data_val = metaAry.data * yscale
        
        pax, = ax.plot(x_axis, data_val, label=metAry['name'])
        
        line_lst.append(pax)
    
    # Set plot limits
    
    # x0, x1
    # xunit, xscale, xlabl
    # ax0
    ax0.set_xlim([x0_scl, x1_scl])
    ax0.set_xlabel(xlabl, fontsize=fontsize)
    ax0.grid(grid)
       
    ax_lst = []
    for ykey in y_lst.keys():
        
        y0 = y_lst[ykey]['y0']
        y1 = y_lst[ykey]['y1']
        ylabl = y_lst[ykey]['ylabl']
        ax = y_lst[ykey]['ax']
        
        ax.set_ylim([y0, y1])
        ax.set_ylabel(ylabl, fontsize=fontsize)
        ax_lst.append(ax)
    
    if legend >= 0:
        ax0.legend(line_lst, [l.get_label() for l in line_lst])
        ax0.legend(loc=legend)
    
    
    # Define the title
    if isinstance(title, str):
        # Use the title string as is
        if fontsize is not None:
            ax0.set_title(title, fontsize=int(fontsize*1.3))
        else:
            ax0.set_title(title)
    elif title:
        t_lst = []
        for metAry in metaArylst:
            t_lst.append(metAry['name'])
        
        title = 'Comparison between ' + ' & '.join(t_lst)
     
    return fig, ax_lst



