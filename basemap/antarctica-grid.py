##
# antactica-grid.py - Created by Timothy Morey on 10/30/2012
#
# This script produces a graphic showing how PISM's horizontal grid is laid
# out and how columns of ice are mapped to processes.
#

from mpl_toolkits.basemap import Basemap #, NetCDFFile

# importing NetCDFFile on the previous line fails, as it is a deprecated object.
# http://matplotlib.1069221.n5.nabble.com/error-installing-basemap-td38910.html
# says that we can import it like so:
from netCDF4 import Dataset as NetCDFFile

import matplotlib.pyplot as plt
import sys
 
nc = NetCDFFile("a10km.cdf1.nc", 'r')
 
# x and y *in the dataset* are only used to determine plotting domain
# dimensions
x = nc.variables['x'][:]
y = nc.variables['y'][:]
width = x.max() - x.min()
height = y.max() - y.min()
 
m = Basemap(width=width,      # width in projection coordinates, in meters
            height=height,      # height
            resolution='l',     # coastline resolution, can be 'l' (low), 'h'
                                # (high) and 'f' (full)
            projection='stere', # stereographic projection
            lat_ts=-71,          # latitude of true scale
            lon_0=90,          # longitude of the plotting domain center
            lat_0=-88)           # latitude of the plotting domain center
 
m.drawlsmask()

width = m.urcrnrx - m.llcrnrx
height = m.urcrnry - m.llcrnry

# draw vertical process boundaries:
for i in range(1,8):
    xcoord = m.llcrnrx + i*width/8
    plt.plot([xcoord, xcoord], [m.llcrnry, m.urcrnry], 'k-')

# draw horizontal process boundaries:
for i in range(1,8):
    ycoord = m.llcrnry + i*height/8
    plt.plot([m.llcrnrx, m.urcrnrx], [ycoord, ycoord], 'k-')

# draw rank labels
for row in range(8):
    for col in range(8):
        plt.text(m.llcrnrx + col*width/8 + width/16,
                 m.llcrnry + row*height/8 + height/16,
                 str(row + col*8),
                 ha='center', va='center')
    
#plt.title(r"Modeled ice surface velocity, m/year")
plt.savefig('antarctica-grid.png')
