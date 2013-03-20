#! /bin/python

##
# build_dataset.py - Created by Timothy Morey on 3/12/2013
#
# This script builds a master dataset from the various cresis grids for
# Greenland or Antarctica.
#


from cresis_grid import CresisGrid
from searise_grid import SeariseGrid

import getopt
import math
import netCDF4
import numpy
import os
import pyproj
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


def createOutfile(outfile, thkgrids, topggrids, fillgrid):
    # Assuming that all grids have the same resolution, and that the thk,topg
    # grids overlap perfectly.

    xmin = sys.float_info.max
    xmax = -sys.float_info.max
    xres = 0.0
    ymin = sys.float_info.max
    ymax = -sys.float_info.max
    yres = 0.0

    cresiscrs = None
    searisecrs = fillgrid.getCRS()

    print thkgrids, len(thkgrids)
    for grid in thkgrids:
        xmin = min(xmin, grid.getMinX())
        xmax = max(xmax, grid.getMaxX())
        xres = grid.cellsize
        ymin = min(ymin, grid.getMinY())
        ymax = max(ymax, grid.getMaxY())
        yres = grid.cellsize
        cresiscrs = grid.getCRS()
    print thkgrids, len(thkgrids)

    searisebounds = [(fillgrid.getMinX(), fillgrid.getMinY()),
                     (fillgrid.getMinX(), fillgrid.getMaxY()),
                     (fillgrid.getMaxX(), fillgrid.getMinY()),
                     (fillgrid.getMaxX(), fillgrid.getMaxY())]
    for pt in searisebounds:
        cx, cy = pyproj.transform(searisecrs, cresiscrs, pt[0], pt[1])
        if cx < xmin:
            xmin -= math.ceil((xmin - cx) / xres) * xres
        if cx > xmax:
            xmax += math.ceil((cx - xmax) / xres) * xres
        if cy < ymin:
            ymin -= math.ceil((ymin - cy) / yres) * yres
        if cy > ymax:
            ymax += math.ceil((cy - ymax) / yres) * yres

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
        print xs, xn, ys, yn
        thk[ys:ys+yn, xs:xs+xn] = data

    for grid in topggrids:
        data = grid.readData()
        xs = int((grid.getMinX() - xmin) / xres)
        xn = int((grid.getMaxX() - grid.getMinX()) / xres) + 1
        ys = int((grid.getMinY() - ymin) / yres)
        yn = int((grid.getMaxY() - grid.getMinY()) / yres) + 1
        print xs, xn, ys, yn
        topg[ys:ys+yn, xs:xs+xn] = data

    root.close()


def fillGaps(datafile, fillfile):
    root = netCDF4.Dataset(outfile, 'a', format='NETCDF3_CLASSIC')
    searise = SeariseGrid(fillfile)
    x = root.variables['x']
    y = root.variables['y']
    thk = root.variables['thk'][:]
    thkfill = -9999.0 # TODO: use _FillValue attribute
    topg = root.variables['topg'][:]
    topgfill = -9999.0 # TODO: use _FillValue attribute

    seariseCRS = searise.getCRS()
    cresisCRS = pyproj.Proj('+proj=stere +lat_ts=70 +lat_0=90 +lon_0=-45 ' +
                            '+k_0=1.0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84')

    for xi in range(len(x)):
        print xi, len(x)
        for yi in range(len(y)):
            if thk.mask[yi,xi] or topg.mask[yi,xi]:
                sx, sy = pyproj.transform(cresisCRS, seariseCRS, x[xi], y[yi])

            if thk.mask[yi,xi]:
                i = searise.interp('thk', sx, sy)
                if i != -9999.0:
                    thk[yi,xi] = i
                    thk.mask[yi,xi] = False

            if topg.mask[yi,xi]:
                i = searise.interp('topg', sx, sy)
                if i != -9999.0:
                    topg[yi,xi] = i
                    topg.mask[yi,xi] = False

        if xi > 0 and xi % 10 == 0:
            print 'saving progress...'
            root.variables['thk'][:] = thk
            root.variables['topg'][:] = topg
            root.close()
            root = netCDF4.Dataset(outfile, 'a', format='NETCDF3_CLASSIC')
            x = root.variables['x']
            y = root.variables['y']
            thk = root.variables['thk'][:]
            topg = root.variables['topg'][:]

    root.variables['thk'][:] = thk
    root.variables['topg'][:] = topg

    root.close()


if __name__ == '__main__':
    indir = '.'
    outfile = './out.nc'
    fillfile = './Greenland1km.nc'
    opts, remainder = getopt.getopt(sys.argv[1:], 'i:o:f:')
    for opt, arg in opts:
        if opt == '-i':
            indir = arg

        elif opt == '-o':
            outfile = arg

        elif opt == '-f':
            fillfile = arg

    # HACK: there is a little overlap between the Jakobshavn and NWCoast data in
    # Greenland, and it looks like Jakobshavn has better data in the region of
    # overlap.  By reversing these lists, we ensure that we put the Jako data
    # into the aggregate file after NWCoast, so we overwrite the NWCoast data.
    thkgrids = [i for i in reversed([CresisGrid(thkfile) for thkfile in findThkFiles(indir)])]
    topggrids = [i for i in reversed([CresisGrid(topgfile) for topgfile in findTopgFiles(indir)])]
    fillgrid = SeariseGrid(fillfile)

    createOutfile(outfile, thkgrids, topggrids, fillgrid)
    fillGaps(outfile, fillfile)
