/**
   large-file.c - Created by Timothy Morey on 10/15/2012

   This file attempts to write out some large CDF-5 netcdf files using Argonne's
   parallel-netcdf library.  We have been unable to do this using pism on ranger,
   but we haven't yet reproduced this on tessy, presumably because we aren't able
   to scale pism large enough.  Hopefully this will allow us to reproduce the
   problem locally.
 */

#include <getopt.h>
#include <math.h>
#include <pnetcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NVARS 4
#define NDIMS 3
#define XNAME "x"
#define XDIM 0
#define XSIZE 1000
#define YNAME "y"
#define YDIM 1
#define YSIZE 1000
#define ZNAME "z"
#define ZDIM 2
#define ZSIZE 1000


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
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  int mpirank, mpisize;

  int ncresult;
  int ncid;
  int varid[NVARS];
  int dimid[NDIMS];

  char c;
  char varname[32];
  char filename[256];

  double* data = 0;
  int i, ix, iy, iz;
  MPI_Offset start[NDIMS], count[NDIMS];

  int localx, localy, localwidth, localheight;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  sprintf(filename, "output.nc");

  while((c = getopt(argc, argv, "o:")) != -1)
  {
    switch(c)
    {
    case 'o':
      strcpy(filename, optarg);  break;
    }
  }

  if(! mpirank) printf("Setting up dimensions...\n");

  ncresult = ncmpi_create(comm, filename, NC_64BIT_DATA, info, &ncid);
  ncerror("ncmpi_create", ncresult);

  ncresult = ncmpi_def_dim(ncid, XNAME, XSIZE, &dimid[XDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_dim(ncid, YNAME, YSIZE, &dimid[YDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_dim(ncid, ZNAME, ZSIZE, &dimid[ZDIM]);  
  ncerror("ncmpi_def_dim", ncresult);

  if(! mpirank) printf("Closing...\n");

  ncresult = ncmpi_close(ncid);
  ncerror("ncmpi_close", ncresult);

  if(! mpirank) printf("Creating some data...\n");
  
  // Distribute our data values in a pism-y way
  GetLocalBounds(XSIZE, YSIZE, mpirank, mpisize,
                   &localx, &localy, &localwidth, &localheight);
  printf("Rank%02d: x=%d, y=%d, width=%d, height=%d\n",
         mpirank, localx, localy, localwidth, localheight);
  start[XDIM] = localx;  start[YDIM] = localy;  start[ZDIM] = 0;
  count[XDIM] = localwidth;  count[YDIM] = localheight;  count[ZDIM] = ZSIZE;
  
  data = (double*)malloc(localwidth * localheight * ZSIZE * sizeof(double));
  for(ix = 0; ix < localwidth; ix++)
    for(iy = 0; iy < localheight; iy++)
      for(iz = 0; iz < ZSIZE; iz++)
        data[iz*localwidth*localheight + iy*localwidth + ix] = 
          iz*XSIZE*YSIZE + (iy+localy)*XSIZE + ix+localx;
  
  for(i = 0; i < NVARS; i++)
  {
    if(! mpirank) printf("Openning...\n");

    ncresult = ncmpi_open(comm, filename, NC_WRITE, info, &ncid);
    ncerror("ncmpi_open", ncresult);

    ncmpi_redef(ncid);

    sprintf(varname, "var%02d", i);
    ncresult = ncmpi_def_var(ncid, varname, NC_DOUBLE, NDIMS, dimid, &varid[i]);
    ncerror("ncmpi_def_var", ncresult);

    ncmpi_enddef(ncid);
    
    if(! mpirank) printf("Writing data...\n");

    ncresult = ncmpi_put_vara_double_all(ncid, varid[i], start, count, data);
    ncerror("ncmpi_put_vara_double", ncresult);
    
    if(! mpirank) printf("Closing file...\n");
    
    ncresult = ncmpi_close(ncid);
    ncerror("ncmpi_close", ncresult);
  }

  if(! mpirank) printf("Done.\n");

  free(data);
  MPI_Finalize();

  return 0;
}
