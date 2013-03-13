
/**
   pncwrite.c - Created by Timothy Morey on 3/4/2013

 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pnetcdf.h>

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

//#define OLD_WRITE_PATTERN 1


int ncerror(const char* function, int result)
{
  if(result)
  {
    fprintf(stderr, "NetCDF function \"%s\" failed with \"%s\"\n",
	    function, ncmpi_strerror(result));
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
  MPI_Offset start[NDIMS], count[NDIMS];  
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
      data2d[ix*localheight + iy] = (ix+localx)*YSIZE + iy+localy;
      for(iz = 0; iz < ZSIZE; iz++)
      {
        z[iz] = iz;
        data3d[ix*localheight*ZSIZE + iy*ZSIZE + iz] = 
          (ix+localx)*YSIZE*ZSIZE + (iy+localy)*ZSIZE + iz;
      }
    }
  }

  MPI_Info_create(&mpiinfo);
  MPI_Info_set(mpiinfo, "use_pism_customizations", "1");

  if(! mpirank) printf("Creating NetCDF file...\n");

  ncresult = ncmpi_create(mpicomm, filename, NC_64BIT_DATA, mpiinfo, &ncid);
  ncerror("ncmpi_create", ncresult);

  if(! mpirank) printf("Setting up dimensions...\n");

  ncresult = ncmpi_def_dim(ncid, TNAME, NC_UNLIMITED, &dimid[TDIM]);
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, TNAME, NC_DOUBLE, 1, &dimid[TDIM], &tvarid);
  ncerror("ncmpi_def_var", ncresult);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing t...\n");
  ncmpi_enddef(ncid);
  start[0] = 0;  count[0] = 1;
  ncresult = ncmpi_put_vara_double_all(ncid, tvarid, start, count, &t);
  ncerror("ncmpi_put_var (t)", ncresult);
  ncmpi_redef(ncid);
#endif

  ncresult = ncmpi_def_dim(ncid, XNAME, XSIZE, &dimid[XDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, XNAME, NC_DOUBLE, 1, &dimid[XDIM], &xvarid);
  ncerror("ncmpi_def_var", ncresult);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing x...\n");
  ncmpi_enddef(ncid);
  start[0] = localx;  count[0] = localwidth;
  ncresult = ncmpi_put_vara_double_all(ncid, xvarid, start, count, x);
  ncerror("ncmpi_put_vara_all (x)", ncresult);
  ncmpi_redef(ncid);
#endif

  ncresult = ncmpi_def_dim(ncid, YNAME, YSIZE, &dimid[YDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, YNAME, NC_DOUBLE, 1, &dimid[YDIM], &yvarid);
  ncerror("ncmpi_def_var", ncresult);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing y...\n");
  ncmpi_enddef(ncid);
  start[0] = localy;  count[0] = localheight;
  ncresult = ncmpi_put_vara_double_all(ncid, yvarid, start, count, y);
  ncerror("ncmpi_put_vara_double_all (y)", ncresult);
  ncmpi_redef(ncid);
#endif

  ncresult = ncmpi_def_dim(ncid, ZNAME, ZSIZE, &dimid[ZDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, ZNAME, NC_DOUBLE, 1, &dimid[ZDIM], &zvarid);
  ncerror("ncmpi_def_var", ncresult);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing z...\n");
  ncmpi_enddef(ncid);
  ncresult = ncmpi_put_var_double_all(ncid, zvarid, z);
  ncerror("ncmpi_put_var_double_all (z)", ncresult);
  ncmpi_redef(ncid);
#endif

  if(! mpirank) printf("Defining variables...\n");

  MPI_Allreduce(&localwidth, &maxwidth, 1, MPI_INT, MPI_MAX, mpicomm);
  MPI_Allreduce(&localheight, &maxheight, 1, MPI_INT, MPI_MAX, mpicomm);
  
  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var3d-%02d", i);

    ncresult = ncmpi_def_var(ncid, varname, NC_DOUBLE, NDIMS, dimid, &var3did[i]);
    ncerror("ncmpi_def_var", ncresult);

    sprintf(varname, "var2d-%02d", i);
    ncresult = ncmpi_def_var(ncid, varname, NC_DOUBLE, NDIMS-1, dimid, &var2did[i]);
  }

  ncmpi_enddef(ncid);

#ifndef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing t...\n");
  start[0] = 0;  count[0] = 1;
  ncresult = ncmpi_put_vara_double_all(ncid, tvarid, start, count, &t);
  ncerror("ncmpi_put_var (t)", ncresult);

  if(! mpirank) printf("Writing x...\n");
  start[0] = localx;  count[0] = localwidth;
  ncresult = ncmpi_put_vara_double_all(ncid, xvarid, start, count, x);
  ncerror("ncmpi_put_vara_all (x)", ncresult);

  if(! mpirank) printf("Writing y...\n");
  start[0] = localy;  count[0] = localheight;
  ncresult = ncmpi_put_vara_double_all(ncid, yvarid, start, count, y);
  ncerror("ncmpi_put_vara_double_all (y)", ncresult);

  if(! mpirank) printf("Writing z...\n");
  ncresult = ncmpi_put_var_double_all(ncid, zvarid, z);
  ncerror("ncmpi_put_var_double_all (z)", ncresult);
#endif

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

    ncresult = ncmpi_put_vara_double_all(ncid, var3did[i], start, count, data3d);
    ncerror("ncmpi_put_vara_double", ncresult);

    if(! mpirank) printf("Writing var2d-%02d...\n", i);

    ncresult = ncmpi_put_vara_double_all(ncid, var2did[i], start, count, data2d);
    ncerror("ncmpi_put_vara_double", ncresult);
  }

  if(! mpirank) printf("Closing file...\n");

  ncmpi_close(ncid);

  if(! mpirank) printf("Done.\n");

  free(data2d);
  free(data3d);
  free(x);
  free(y);
  free(z);

  MPI_Finalize();

  return 0;
}
