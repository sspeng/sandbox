#! /bin/python

##
# build_dataset.py - Created by Timothy Morey on 3/12/2013
#
# This script builds a master dataset from the various cresis grids for
# Greenland or Antarctica.
#


from cresis_grid import CresisGrid

import getopt
import netCDF4
import numpy
import os
import sys


def findThkFiles(indir):
    thkfiles = []
    for f in os.listdir(indir):
        path = os.path.join(indir, f)
        if os.path.isdir(path):
            path = os.path.join(path, 'grids')
            if os.path.isdir(path):
                for gridfile in os.listdir(path):
                    if 'thickness' in gridfile.lower() and '.txt' in gridfile:
                        thkfiles.append(os.path.join(path, gridfile))

    return thkfiles


def findTopgFiles(indir):
    topgfiles = []
    for f in os.listdir(indir):
        path = os.path.join(indir, f)
        if os.path.isdir(path):
            path = os.path.join(path, 'grids')
            if os.path.isdir(path):
                for gridfile in os.listdir(path):
                    if 'bottom' in gridfile.lower() and '.txt' in gridfile:
                        topgfiles.append(os.path.join(path, gridfile))

    return topgfiles


def createOutfile(outfile, thkgrids, topggrids):
    # Assuming that all grids have the same resolution, and that the thk,topg
    # grids overlap perfectly.

    xmin = sys.float_info.max
    xmax = -sys.float_info.max
    xres = 0.0
    ymin = sys.float_info.max
    ymax = -sys.float_info.max
    yres = 0.0

    for grid in thkgrids:
        xmin = min(xmin, grid.getMinX())
        xmax = max(xmax, grid.getMaxX())
        xres = grid.cellsize
        ymin = min(ymin, grid.getMinY())
        ymax = max(ymax, grid.getMaxY())
        yres = grid.cellsize

    root = netCDF4.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')

    # Define dimensions
    xdim = root.createDimension('x', int((xmax - xmin) / xres) + 1)
    ydim = root.createDimension('y', int((ymax - ymin) / yres) + 1)

    # Define variables
    xvar = root.createVariable('x', 'f4', ('x',))
    yvar = root.createVariable('y', 'f4', ('y',))
    topg = root.createVariable('topg', 'f4', ('y', 'x'), fill_value=-9999.0)
    thk = root.createVariable('thk', 'f4', ('y', 'x'), fill_value=-9999.0)

    # Write dimension coordinates
    xvar[:] = numpy.arange(xmin, xmax + xres, xres)
    yvar[:] = numpy.arange(ymin, ymax + yres, yres)

    for grid in thkgrids:
        data = grid.readData()
        xs = int((grid.getMinX() - xmin) / xres)
        xn = int((grid.getMaxX() - grid.getMinX()) / xres) + 1
        ys = int((grid.getMinY() - ymin) / yres)
        yn = int((grid.getMaxY() - grid.getMinY()) / yres) + 1
        thk[ys:ys+yn, xs:xs+xn] = data

    for grid in topggrids:
        data = grid.readData()
        xs = int((grid.getMinX() - xmin) / xres)
        xn = int((grid.getMaxX() - grid.getMinX()) / xres) + 1
        ys = int((grid.getMinY() - ymin) / yres)
        yn = int((grid.getMaxY() - grid.getMinY()) / yres) + 1
        topg[ys:ys+yn, xs:xs+xn] = data

    root.close()


def addTopgGridToOutput(outfile, cresisGrid):
    root = netCDF4.Dataset(outfile, 'a', format='NETCDF3_CLASSIC')
    topg = root.variables['topg']
    data = cresisGrid.readData()

if __name__ == '__main__':
    indir = '.'
    outfile = './out.nc'
    opts, remainder = getopt.getopt(sys.argv[1:], 'i:o:')
    for opt, arg in opts:
        if opt == '-i':
            indir = arg
            pass
        elif opt == '-o':
            outfile = arg
            pass

    # HACK: there is a little overlap between the Jakobshavn and NWCoast data in
    # Greenland, and it looks like Jakobshavn has better data in the region of
    # overlap.  By reversing these lists, we ensure that we put the Jako data
    # into the aggregate file after NWCoast, so we overwrite the NWCoast data.
    thkgrids = reversed([CresisGrid(thkfile) for thkfile in findThkFiles(indir)])
    topggrids = reversed([CresisGrid(topgfile) for topgfile in findTopgFiles(indir)])

    createOutfile(outfile, thkgrids, topggrids)
