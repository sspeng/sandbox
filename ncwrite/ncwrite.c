/**
   ncwrite.c - Created by Timothy Morey on 10/26/2012

   This file attempts to reproduce the write error we have been seeing in PISM
   without all of the confusion associated with PISM.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf_par.h>
#include <netcdf.h>

#define NVARS 3
#define NDIMS 4
#define XNAME "x"
#define XDIM 1
#define XSIZE 1501
#define YNAME "y"
#define YDIM 2
#define YSIZE 2801
#define ZNAME "z"
#define ZDIM 3
#define ZSIZE 401
#define TNAME "t"
#define TDIM 0


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

int GetLocalBounds(int gwidth, int gheight, int rank, int size, 
		   int* x, int* y, int* width, int* height)
{
  int retval = 0;
  int sizex = 0;
  int sizey = 0;

  // Using PISM's algorithm for distributing gwidth columns and gheight rows of
  // grid cells over size processors:
  sizex = (int)(0.5 + sqrt(((double)gwidth)*((double)size)/((double)gheight)));
  if(sizex == 0)
    sizex = 1;

  while(sizex > 0)
  {
    sizey = size / sizex;
    if(size == sizex*sizey)
      break;
    sizex--;
  }

  if(gwidth > gheight && sizex < sizey)
  {
    int temp = sizex;
    sizex = sizey;
    sizey = temp;
  }

  if(sizex == 0 || 
     sizey == 0 ||
     (gwidth / sizex) < 1 ||
     (gheight / sizey) < 1)
  {
    printf("Failed to distribute %d columns and %d rows over %d processors\n", 
	   gwidth, gheight, size);
    retval = 1;
  }
  else
  {
    int extrac = gwidth % sizex;
    int extrar = gheight % sizey;
    int r = rank / sizex;
    int c = rank % sizex;

    *width = gwidth / sizex;
    if(c < extrac)
      *width += 1;

    *height = gheight / sizey;
    if(r < extrar)
      *height += 1;

    *x = (*width) * c + (c < extrac ? c : extrac);
    *y = (*height) * r + (r < extrar ? r : extrar);
  }
  
  return retval;
}

int main(int argc, char* argv[])
{
  char c;
  int ix, iy, iz, i;
  MPI_Comm mpicomm;
  MPI_Info mpiinfo;
  int mpirank;
  int mpisize;
  double *data3d, *data2d, *x, *y, *z, t;
  int localx, localy, localwidth, localheight;
  int maxwidth, maxheight;
  int ncresult;
  int ncid, tvarid, xvarid, yvarid, zvarid;
  int dimid[NDIMS], var3did[NVARS], var2did[NVARS];
  size_t start[NDIMS], count[NDIMS], chunksize[NDIMS];  
  const char* filename = "output.nc";
  char varname[32];
  
  mpicomm = MPI_COMM_WORLD;
  mpiinfo = MPI_INFO_NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(mpicomm, &mpisize);
  MPI_Comm_rank(mpicomm, &mpirank);

  if(! mpirank) printf("Creating some data...\n");
  
  // Distribute our data values in a pism-y way
  GetLocalBounds(XSIZE, YSIZE, mpirank, mpisize,
                 &localx, &localy, &localwidth, &localheight);
  printf("Rank%02d: x=%d, y=%d, width=%d, height=%d\n",
         mpirank, localx, localy, localwidth, localheight);
  
  data2d = (double*)malloc(localwidth * localheight * sizeof(double));
  data3d = (double*)malloc(localwidth * localheight * ZSIZE * sizeof(double));
  x = (double*)malloc(localwidth * sizeof(double));
  y = (double*)malloc(localheight * sizeof(double));
  z = (double*)malloc(ZSIZE * sizeof(double));
  t = 0.0;
  for(ix = 0; ix < localwidth; ix++)
  {
    x[ix] = ix + localx;
    for(iy = 0; iy < localheight; iy++)
    {
      y[iy] = iy + localy;
      data2d[ix*localheight + iy] = (ix+localx)*localheight + iy+localy;
      for(iz = 0; iz < ZSIZE; iz++)
      {
        z[iz] = iz;
        data3d[ix*localheight*ZSIZE + iy*ZSIZE + iz] = 
          (ix+localx)*YSIZE*ZSIZE + (iy+localy)*ZSIZE + iz;
      }
    }
  }

  if(! mpirank) printf("Creating NetCDF file...\n");

  ncresult = nc_create_par(filename, NC_MPIIO | NC_NETCDF4, 
                           mpicomm, mpiinfo, &ncid);
  ncerror("nc_create_par", ncresult);

  if(! mpirank) printf("Setting up dimensions...\n");

  ncresult = nc_def_dim(ncid, TNAME, NC_UNLIMITED, &dimid[TDIM]);
  ncerror("nc_def_dim", ncresult);

  ncresult = nc_def_var(ncid, TNAME, NC_DOUBLE, 1, &dimid[TDIM], &tvarid);
  ncerror("nc_def_var", ncresult);

  if(! mpirank) printf("Writing t...\n");

  nc_enddef(ncid);

  ncresult = nc_var_par_access(ncid, tvarid, NC_COLLECTIVE);
  ncerror("nc_var_par_access", ncresult);

  start[0] = 0;  count[0] = 1;
  ncresult = nc_put_vara_double(ncid, tvarid, start, count, &t);
  ncerror("nc_put_var (t)", ncresult);

  nc_redef(ncid);

  ncresult = nc_def_dim(ncid, XNAME, XSIZE, &dimid[XDIM]);  
  ncerror("nc_def_dim", ncresult);

  ncresult = nc_def_var(ncid, XNAME, NC_DOUBLE, 1, &dimid[XDIM], &xvarid);
  ncerror("nc_def_var", ncresult);

  if(! mpirank) printf("Writing x...\n");

  nc_enddef(ncid);

  ncresult = nc_var_par_access(ncid, xvarid, NC_COLLECTIVE);
  ncerror("nc_var_par_access", ncresult);

  start[0] = localx;  count[0] = localwidth;
  ncresult = nc_put_vara_double(ncid, xvarid, start, count, x);
  ncerror("nc_put_vara_all (x)", ncresult);

  nc_redef(ncid);

  ncresult = nc_def_dim(ncid, YNAME, YSIZE, &dimid[YDIM]);  
  ncerror("nc_def_dim", ncresult);

  ncresult = nc_def_var(ncid, YNAME, NC_DOUBLE, 1, &dimid[YDIM], &yvarid);
  ncerror("nc_def_var", ncresult);

  if(! mpirank) printf("Writing y...\n");

  nc_enddef(ncid);

  ncresult = nc_var_par_access(ncid, yvarid, NC_COLLECTIVE);
  ncerror("nc_var_par_access", ncresult);

  start[0] = localy;  count[0] = localheight;
  ncresult = nc_put_vara_double(ncid, yvarid, start, count, y);
  ncerror("nc_put_vara_all (y)", ncresult);

  nc_redef(ncid);

  ncresult = nc_def_dim(ncid, ZNAME, ZSIZE, &dimid[ZDIM]);  
  ncerror("nc_def_dim", ncresult);

  ncresult = nc_def_var(ncid, ZNAME, NC_DOUBLE, 1, &dimid[ZDIM], &zvarid);
  ncerror("nc_def_var", ncresult);

  if(! mpirank) printf("Writing z...\n");

  nc_enddef(ncid);

  ncresult = nc_var_par_access(ncid, zvarid, NC_COLLECTIVE);
  ncerror("nc_var_par_access", ncresult);

  ncresult = nc_put_var_double(ncid, zvarid, z);
  ncerror("nc_put_var (z)", ncresult);

  nc_redef(ncid);

  if(! mpirank) printf("Defining variables...\n");

  MPI_Allreduce(&localwidth, &maxwidth, 1, MPI_INT, MPI_MAX, mpicomm);
  MPI_Allreduce(&localheight, &maxheight, 1, MPI_INT, MPI_MAX, mpicomm);
  
  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var3d-%02d", i);
    ncresult = nc_def_var(ncid, varname, NC_DOUBLE, NDIMS, dimid, &var3did[i]);
    ncerror("nc_def_var", ncresult);

    chunksize[TDIM] = 1;
    chunksize[XDIM] = maxwidth;
    chunksize[YDIM] = maxheight;
    chunksize[ZDIM] = ZSIZE;
    ncresult = nc_def_var_chunking(ncid, var3did[i], NC_CHUNKED, chunksize);
    ncerror("nc_def_var_chunking", ncresult);

    sprintf(varname, "var2d-%02d", i);
    ncresult = nc_def_var(ncid, varname, NC_DOUBLE, NDIMS-1, dimid, &var2did[i]);

    ncresult = nc_def_var_chunking(ncid, var2did[i], NC_CHUNKED, chunksize);
    ncerror("nc_def_var_chunking", ncresult);
  }

  nc_enddef(ncid);

  if(! mpirank) printf("Writing variable data...\n");

  for(i = 0; i < NVARS; i++)
  {
    start[XDIM] = localx; 
    start[YDIM] = localy;  
    start[ZDIM] = 0;  
    start[TDIM] = 0;
    
    count[XDIM] = localwidth;  
    count[YDIM] = localheight;  
    count[ZDIM] = ZSIZE;
    count[TDIM] = 1;

    if(! mpirank) printf("Writing var3d-%02d...\n", i);

    ncresult = nc_var_par_access(ncid, var3did[i], NC_COLLECTIVE);
    ncerror("nc_var_par_access", ncresult);

    ncresult = nc_put_vara_double(ncid, var3did[i], start, count, data3d);
    ncerror("nc_put_vara_double", ncresult);

    if(! mpirank) printf("Writing var2d-%02d...\n", i);

    ncresult = nc_var_par_access(ncid, var2did[i], NC_COLLECTIVE);
    ncerror("nc_var_par_access", ncresult);

    ncresult = nc_put_vara_double(ncid, var2did[i], start, count, data2d);
    ncerror("nc_put_vara_double", ncresult);
  }

  if(! mpirank) printf("Closing file...\n");

  nc_close(ncid);

  if(! mpirank) printf("Done.\n");

  free(data2d);
  free(data3d);
  free(x);
  free(y);
  free(z);

  MPI_Finalize();

  return 0;
}
