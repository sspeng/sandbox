/**
   pism-latlon.c - Created by Timothy Morey on 10/11/2012

   This file provides a simple tool to preprocess netcdf files that might be 
   used as bootstrap-input for PISM.  It will calculate the lat/lon variables
   from the x/y values, and store them in the file.
 */

#include <netcdf.h>
#include <proj_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define XDIM "x"
#define YDIM "y"

#define XVAR "x"
#define YVAR "y"
#define LATVAR "lat"
#define LONVAR "lon"

#define GRID_MAPPING "mapping"
#define LAT_UNITS "degreeN"
#define LON_UNITS "degreeE"
#define LAT_LONG_NAME "Latitude"
#define LAT_STD_NAME "latitude"
#define LON_LONG_NAME "Longitude"
#define LON_STD_NAME "longitude"

int ncerror(const char* function, int result)
{
  if(result)
  {
    fprintf(stderr, "NetCDF function \"%s\" failed with \"%s\"\n",
	    function, nc_strerror(result));
    exit(1);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  int ncid, xvarid, yvarid, latvarid, lonvarid;
  int xydims[2];
  size_t xlen, ylen;
  int res;
  double *x, *y, *lat, *lon;
  int ix, iy;
  double tempx, tempy;

  projPJ pj_pism, pj_wgs84;

  res = nc_open(argv[1], NC_WRITE, &ncid); ncerror("nc_open", res);
  res = nc_inq_dimid(ncid, XDIM, &xydims[0]); ncerror("nc_inq_dimid", res);
  res = nc_inq_dimlen(ncid, xydims[0], &xlen); ncerror("nc_inq_dimlen", res);
  res = nc_inq_dimid(ncid, YDIM, &xydims[1]); ncerror("nc_inq_dimid", res);
  res = nc_inq_dimlen(ncid, xydims[1], &ylen); ncerror("nc_inq_dimlen", res);
  res = nc_inq_varid(ncid, XVAR, &xvarid); ncerror("nc_inq_varid", res);
  res = nc_inq_varid(ncid, YVAR, &yvarid); ncerror("nc_inq_varid", res);

  if(res = nc_inq_varid(ncid, LATVAR, &latvarid))
  {
    printf("Defining lat...\n");
    res = nc_redef(ncid);

    res = nc_def_var(ncid, LATVAR, NC_DOUBLE, 2, xydims, &latvarid);
    ncerror("nc_def_var", res);

    res = nc_put_att_text(ncid, latvarid, "grid_mapping", 
			  strlen(GRID_MAPPING), GRID_MAPPING);
    res = nc_put_att_text(ncid, latvarid, "long_name",
			  strlen(LAT_LONG_NAME), LAT_LONG_NAME);
    res = nc_put_att_text(ncid, latvarid, "standard_name",
			  strlen(LAT_STD_NAME), LAT_STD_NAME);
    res = nc_put_att_text(ncid, latvarid, "units",
			  strlen(LAT_UNITS), LAT_UNITS);
  }

  if(res = nc_inq_varid(ncid, LONVAR, &lonvarid))
  {
    printf("Defining lon...\n");
    res = nc_redef(ncid);

    res = nc_def_var(ncid, LONVAR, NC_DOUBLE, 2, xydims, &lonvarid);
    ncerror("nc_def_var", res);

    res = nc_put_att_text(ncid, lonvarid, "grid_mapping", 
			  strlen(GRID_MAPPING), GRID_MAPPING);
    res = nc_put_att_text(ncid, lonvarid, "long_name",
			  strlen(LON_LONG_NAME), LON_LONG_NAME);
    res = nc_put_att_text(ncid, lonvarid, "standard_name",
			  strlen(LON_STD_NAME), LON_STD_NAME);
    res = nc_put_att_text(ncid, lonvarid, "units",
			  strlen(LON_UNITS), LON_UNITS);
  }

  res = nc_enddef(ncid);

  x = (double*)malloc(xlen * sizeof(double));
  y = (double*)malloc(ylen * sizeof(double));
  lat = (double*)malloc(xlen * ylen * sizeof(double));
  lon = (double*)malloc(xlen * ylen * sizeof(double));

  printf("Retrieving x values...\n");
  res = nc_get_var_double(ncid, xvarid, x); ncerror("nc_get_var_double", res);

  printf("Retrieving y values...\n");
  res = nc_get_var_double(ncid, yvarid, y); ncerror("nc_get_var_double", res);

  printf("Initializing proj4...\n");
  pj_pism = pj_init_plus("+proj=stere +lat_ts=71 +lat_0=90 +lon_0=-39 +k_0=1.0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84");
  pj_wgs84 = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");

  printf("Converting...\n");
  for(ix = 0; ix < xlen; ix++)
  {
    for(iy = 0; iy < ylen; iy++)
    {
      tempx = x[ix];
      tempy = y[iy];
      pj_transform(pj_pism, pj_wgs84, 1, 1, &tempx, &tempy, NULL);
      lon[iy * xlen + ix] = tempx / DEG_TO_RAD;
      lat[iy * xlen + ix] = tempy / DEG_TO_RAD;
    }
  }

  printf("Storing lat/lon...\n");
  res = nc_put_var_double(ncid, lonvarid, lon); ncerror("nc_put_var_double", res);
  res = nc_put_var_double(ncid, latvarid, lat); ncerror("nc_put_var_double", res);

  printf("Closing...\n");
  res = nc_close(ncid); ncerror("nc_close", res);

  printf("Done.\n");

  return 0;
}
