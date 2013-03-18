#! /bin/python

##
# searise_grid.py - Created by Timothy Morey on 3/16/2013
#
# This file contains a definition for the SeariseGrid class, which represents a
# SeaRISE dataset stored as a NetCDF file.
#


import netCDF4
import numpy
import pyproj
import sys


class SeariseGrid:

    def __init__(self, filename):
        self.filename = filename
        self.cachedvars = dict()
        self.readMetadata()

    def getCRS(self):
        # TODO: parse the CRS out of the metadata
        return pyproj.Proj('+proj=stere +lat_ts=71 +lat_0=90 +lon_0=-39 ' +
                           '+k_0=1.0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84')

    def getMinX(self):
        return self.minx

    def getMaxX(self):
        return self.maxx

    def getMinY(self):
        return self.miny

    def getMaxY(self):
        return self.maxy

    def interp(self, varname, x, y):
        if not 'x' in self.cachedvars:
            self.cachedvars['x'] = self.readVar('x')
        if not 'y' in self.cachedvars:
            self.cachedvars['y'] = self.readVar('y')
        if not varname in self.cachedvars:
            self.cachedvars[varname] = self.readVar(varname)

        xvar = self.cachedvars['x']
        yvar = self.cachedvars['y']
        var = self.cachedvars[varname]

        xi1 = None
        xi2 = None
        yi1 = None
        yi2 = None
        for i in range(1, len(xvar)):
            if xvar[i-1] <= x and x <= xvar[i]:
                xi1, xi2 = i-1, i
                x1, x2 = xvar[xi1], xvar[xi2]

        for i in range(1, len(yvar)):
            if yvar[i-1] <= y and y <= yvar[i]:
                yi1, yi2 = i-1, i
                y1, y2 = yvar[yi1], yvar[yi2]

        value = -9999.0
        if (xi1 is not None and
            xi2 is not None and
            yi1 is not None and
            yi2 is not None):
            #print xi1, xi2, yi1, yi2, len(var)
            values = [var[0, yi1, xi1],
                      var[0, yi1, xi2],
                      var[0, yi2, xi2],
                      var[0, yi2, xi1]]
            weights = [((x2-x)*(y2-y))/((x2-x1)*(y2-y1)),
                       ((x2-x)*(y-y1))/((x2-x1)*(y2-y1)),
                       ((x-x1)*(y-y1))/((x2-x1)*(y2-y1)),
                       ((x-x1)*(y2-y))/((x2-x1)*(y2-y1))]
            value = 0.0
            for i in range(4):
                value += values[i] * weights[i]

        #print 'interp', x, y, value

        return value

    def readMetadata(self):
        nc = netCDF4.Dataset(self.filename, 'r')

        self.minx = sys.float_info.max
        self.maxx = -sys.float_info.max
        self.miny = sys.float_info.max
        self.maxy = -sys.float_info.max

        xvar = nc.variables['x']
        self.nx = len(xvar)
        for x in xvar:
            self.minx = min(x, self.minx)
            self.maxx = max(x, self.maxx)

        yvar = nc.variables['y']
        self.ny = len(yvar)
        for y in yvar:
            self.miny = min(y, self.miny)
            self.maxy = max(y, self.maxy)

        nc.close()

    def readVar(self, varname):
        nc = netCDF4.Dataset(self.filename, 'r')
        vara = nc.variables[varname][:]
        nc.close()
        return vara
