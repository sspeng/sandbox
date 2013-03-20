#! /bin/python

##
# cresis_grid.py - Created by Timothy Morey on 3/12/2013
#
# This file contains a definition for the CresisGrid class, which represents a
# single 2d grid of CReSIS data.
#


import numpy
import pyproj


class CresisGrid:

    def __init__(self, filename):
        self.filename = filename
        self.readMetadata()

    def getCRS(self):
        # TODO: parse the crs out of a prj file
        return pyproj.Proj('+proj=stere +lat_ts=70 +lat_0=90 +lon_0=-45 ' +
                           '+k_0=1.0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84')


    def getMinX(self):
        return self.xllcorner

    def getMaxX(self):
        return self.xllcorner + self.cellsize*(self.ncols - 1)

    def getMinY(self):
        return self.yllcorner

    def getMaxY(self):
        return self.yllcorner + self.cellsize*(self.nrows - 1)

    def readMetadata(self):
        self.ncols = None
        self.nrows = None
        self.xllcorner = None
        self.yllcorner = None
        self.cellsize = None
        self.nodata = None

        f = open(self.filename, 'r')

        tokens = f.readline().split()
        if 'ncols' in tokens[0]:
            self.ncols = int(tokens[1])

        tokens = f.readline().split()
        if 'nrows' in tokens[0]:
            self.nrows = int(tokens[1])

        tokens = f.readline().split()
        if 'xllcorner' in tokens[0]:
            self.xllcorner = float(tokens[1])

        tokens = f.readline().split()
        if 'yllcorner' in tokens[0]:
            self.yllcorner = float(tokens[1])

        tokens = f.readline().split()
        if 'cellsize' in tokens[0]:
            self.cellsize = float(tokens[1])

        tokens = f.readline().split()
        if 'NODATA_value' in tokens[0]:
            self.nodata = int(tokens[1])

        f.close()

    def readData(self):
        data = numpy.empty((self.nrows, self.ncols))
        f = open(self.filename, 'r')
        
        # skip header lines
        for i in range(6):
            f.readline()

        row = self.nrows - 1
        while row >= 0:
            line = f.readline().split()
            col = 0
            for val in line:
                data[row, col] = float(val)
                col += 1
            row -= 1

        f.close()
        return data

    def __str__(self):
        return unicode(self).encode('utf-8')

    def __unicode__(self):
        s = self.filename + ': '
        s += 'ncols=' + str(self.ncols) + ', '
        s += 'nrows=' + str(self.nrows) + ', '
        s += 'xllcorner=' + str(self.xllcorner) + ', '
        s += 'yllcorner=' + str(self.yllcorner) + ', '
        s += 'cellsize=' + str(self.cellsize) + ', '
        s += 'nodata_value=' + str(self.nodata)
        return s
