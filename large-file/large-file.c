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
#define NDIMS 4
#define XNAME "x"
#define XDIM 1
#define YNAME "y"
#define YDIM 2
#define ZNAME "z"
#define ZDIM 3
#define TNAME "t"
#define TDIM 0


typedef enum
{
  Pismy     = 0x1,
  Alternate = 0x2,
  Both      = 0x3
} WriteMethod;


MPI_Comm g_mpicomm;
MPI_Info g_mpiinfo;
int g_mpirank;
int g_mpisize;
double *g_data, *g_x, *g_y, *g_z, g_t;
int g_localx, g_localy, g_localwidth, g_localheight;
int g_xsize = 1501;
int g_ysize = 2801;
int g_zsize = 401;


int DoAlternateWrite(const char* filename);
int DoPismyWrite(const char* filename);
int GetLocalBounds(int gwidth, int gheight, int rank, int size, 
		   int* x, int* y, int* width, int* height);
int NCError(const char* function, int result);
void PrintInfo(MPI_Info info, const char* label);


int main(int argc, char* argv[])
{
  char c;
  int ix, iy, iz, i;
  char basename[256], filename[256];
  int numtries = 1;
  WriteMethod method = Both;
  int pause = 0;
  
  g_mpicomm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(g_mpicomm, &g_mpisize);
  MPI_Comm_rank(g_mpicomm, &g_mpirank);

  MPI_Info_create(&g_mpiinfo);
  // TODO: set some hints...

  strcpy(basename, "output.nc");

  while((c = getopt(argc, argv, "i:o:m:p:x:y:z:")) != -1)
  {
    switch(c)
    {
    case 'i':
      // Number of test iterations
      numtries = atoi(optarg);
      break;

    case 'o':
      // Output file base name
      strcpy(basename, optarg);
      break;

    case 'm':
      // Write method
      if(0 == strcmp("pimsy", optarg) || 
         0 == strcmp("pism", optarg) ||
         0 == strcmp("1", optarg))
      {
        method = Pismy;
      }
      else if(0 == strcmp("alternate", optarg) ||
              0 == strcmp("alt", optarg) ||
              0 == strcmp("2", optarg))
      {
        method = Alternate;
      }
      else
      {
        method = Both;
      }
      break;

    case 'p':
      // Pause between iterations, in seconds
      pause = atoi(optarg);
      break;

    case 'x':
      // Number of data points in the x dimension
      g_xsize = atoi(optarg);
      break;

    case 'y':
      // Number of data points in the y dimension
      g_ysize = atoi(optarg);
      break;

    case 'z':
      // Number of data points in the z dimension
      g_zsize = atoi(optarg);
      break;
    }
  }

  if(! g_mpirank)
  {
    printf("Using a %dx%dx%d grid.\n", g_xsize, g_ysize, g_zsize);
    printf("Running %d test iterations.\n", numtries);
    printf("Using \"%s\" as the base filename.\n", basename);
    printf("Pausing for %d s between tests.\n", pause);
    if(method == Both)
      printf("Testing both write methods.\n");
    else if(method == Pismy)
      printf("Testing pismy write method..\n");
    else if(method == Alternate)
      printf("Testing alternate write method.\n");
  }

  if(! g_mpirank) printf("Creating some data...\n");
  
  // Distribute our data values in a pism-y way
  GetLocalBounds(g_xsize, g_ysize, g_mpirank, g_mpisize,
                 &g_localx, &g_localy, &g_localwidth, &g_localheight);
  printf("Rank%02d: x=%d, y=%d, width=%d, height=%d\n",
         g_mpirank, g_localx, g_localy, g_localwidth, g_localheight);
  
  g_data = (double*)malloc(g_localwidth * g_localheight * g_zsize * sizeof(double));
  g_x = (double*)malloc(g_localwidth * sizeof(double));
  g_y = (double*)malloc(g_localheight * sizeof(double));
  g_z = (double*)malloc(g_zsize * sizeof(double));
  g_t = 0.0;
  for(ix = 0; ix < g_localwidth; ix++)
  {
    g_x[ix] = ix;
    for(iy = 0; iy < g_localheight; iy++)
    {
      g_y[iy] = iy;
      for(iz = 0; iz < g_zsize; iz++)
      {
        g_z[iz] = iz;
        g_data[iz*g_localwidth*g_localheight + iy*g_localwidth + ix] = 
          iz*g_xsize*g_ysize + (iy+g_localy)*g_xsize + ix+g_localx;
      }
    }
  }

  for(i = 0; i < numtries; i++)
  {
    if(method & Pismy)
    {
      sprintf(filename, "%d.pismy.%s", i, basename);
      DoPismyWrite(filename);

      if(i < numtries - 1 ||
         method & Alternate)
        sleep(pause);
    }

    if(method & Alternate)
    {
      sprintf(filename, "%d.alt.%s", i, basename);
      DoAlternateWrite(filename);

      if(i < numtries - 1)
        sleep(pause);
    }
  }

  free(g_data);
  MPI_Info_free(&g_mpiinfo);
  MPI_Finalize();

  return 0;
}

int NCError(const char* function, int result)
{
  if(result)
  {
    fprintf(stderr, "NetCDF function \"%s\" failed with \"%s\"\n",
	    function, ncmpi_strerror(result));
    exit(1);
  }

  return 0;
}

int DoAlternateWrite(const char* filename)
{
  int ncresult;
  int ncid;
  int varid[NVARS], xvarid, yvarid, zvarid, tvarid;
  int dimid[NDIMS];
  MPI_Offset start[NDIMS], count[NDIMS];
  double starttime, elapsed;
  char varname[32];
  int i;

  if(! g_mpirank) printf("Beginning tests of alternate write\n");
  starttime = MPI_Wtime();

  if(! g_mpirank) printf("Creating NetCDF file...\n");

  if(! g_mpirank) PrintInfo(g_mpiinfo, "Requested hints");

  ncresult = ncmpi_create(g_mpicomm, filename, NC_64BIT_DATA, g_mpiinfo, &ncid);
  NCError("ncmpi_create", ncresult);

  if(! g_mpirank)
  {
    MPI_Info active;
    ncmpi_get_file_info(ncid, &active);
    PrintInfo(active, "Active hints");
  }

  if(! g_mpirank) printf("Setting up dimensions...\n");

  ncresult = ncmpi_def_dim(ncid, TNAME, NC_UNLIMITED, &dimid[TDIM]);
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_dim(ncid, XNAME, g_xsize, &dimid[XDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_dim(ncid, YNAME, g_ysize, &dimid[YDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_dim(ncid, ZNAME, g_zsize, &dimid[ZDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  if(! g_mpirank) printf("Defining variables...\n");

  ncresult = ncmpi_def_var(ncid, TNAME, NC_DOUBLE, 1, &dimid[TDIM], &tvarid);
  NCError("ncmpi_def_var", ncresult);

  ncresult = ncmpi_def_var(ncid, XNAME, NC_DOUBLE, 1, &dimid[XDIM], &xvarid);
  NCError("ncmpi_def_var", ncresult);

  ncresult = ncmpi_def_var(ncid, YNAME, NC_DOUBLE, 1, &dimid[YDIM], &yvarid);
  NCError("ncmpi_def_var", ncresult);

  ncresult = ncmpi_def_var(ncid, ZNAME, NC_DOUBLE, 1, &dimid[ZDIM], &zvarid);
  NCError("ncmpi_def_var", ncresult);

  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var%02d", i);
    ncresult = ncmpi_def_var(ncid, varname, NC_DOUBLE, NDIMS, dimid, &varid[i]);
    NCError("ncmpi_def_var", ncresult);
  }

  ncmpi_enddef(ncid);

  if(! g_mpirank) printf("Writing data...\n");
  
  if(! g_mpirank) printf("Writing x...\n");
  ncmpi_enddef(ncid);
  start[0] = g_localx;  count[0] = g_localwidth;
  ncresult = ncmpi_put_vara_double_all(ncid, xvarid, start, count, g_x);
  NCError("ncmpi_put_vara_double_all (x)", ncresult);

  if(! g_mpirank) printf("Writing y...\n");
  ncmpi_enddef(ncid);
  start[0] = g_localy;  count[0] = g_localheight;
  ncresult = ncmpi_put_vara_double_all(ncid, yvarid, start, count, g_y);
  NCError("ncmpi_put_vara_double_all (y)", ncresult);
  
  if(! g_mpirank) printf("Writing z...\n");
  ncresult = ncmpi_put_var_double_all(ncid, zvarid, g_z);
  NCError("ncmpi_put_var_double_all (z)", ncresult);

  if(! g_mpirank) printf("Writing t...\n");
  start[0] = 0;  count[0] = 1;
  ncresult = ncmpi_put_vara_double_all(ncid, tvarid, start, count, &g_t);
  NCError("ncmpi_put_vara_double_all (t)", ncresult);

  for(i = 0; i < NVARS; i++)
  {
    if(! g_mpirank) printf("Writing var%02d...\n", i);

    start[XDIM] = g_localx; 
    start[YDIM] = g_localy;  
    start[ZDIM] = 0;  
    start[TDIM] = 0;
    
    count[XDIM] = g_localwidth;  
    count[YDIM] = g_localheight;  
    count[ZDIM] = g_zsize;
    count[TDIM] = 1;

    ncresult = ncmpi_put_vara_double_all(ncid, varid[i], start, count, g_data);
    NCError("ncmpi_put_vara_double_all", ncresult);
  }

  elapsed = MPI_Wtime() - starttime;

  if(! g_mpirank) printf("Done.  %f s\n", elapsed);

  return 0;
}

int DoPismyWrite(const char* filename)
{
  int ncresult;
  int ncid;
  int varid[NVARS], xvarid, yvarid, zvarid, tvarid;
  int dimid[NDIMS];
  MPI_Offset start[NDIMS], count[NDIMS];
  double starttime, elapsed;
  char varname[32];
  int i;

  if(! g_mpirank) printf("Beginning tests of pism-ish write\n");
  starttime = MPI_Wtime();

  if(! g_mpirank) printf("Creating NetCDF file...\n");

  ncresult = ncmpi_create(g_mpicomm, filename, NC_64BIT_DATA, g_mpiinfo, &ncid);
  NCError("ncmpi_create", ncresult);

  if(! g_mpirank) printf("Setting up dimensions...\n");

  ncresult = ncmpi_def_dim(ncid, TNAME, NC_UNLIMITED, &dimid[TDIM]);
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, TNAME, NC_DOUBLE, 1, &dimid[TDIM], &tvarid);
  NCError("ncmpi_def_var", ncresult);

  if(! g_mpirank) printf("Writing t...\n");
  ncmpi_enddef(ncid);
  start[0] = 0;  count[0] = 1;
  ncresult = ncmpi_put_vara_double_all(ncid, tvarid, start, count, &g_t);
  NCError("ncmpi_put_var (t)", ncresult);
  ncmpi_redef(ncid);

  ncresult = ncmpi_def_dim(ncid, XNAME, g_xsize, &dimid[XDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, XNAME, NC_DOUBLE, 1, &dimid[XDIM], &xvarid);
  NCError("ncmpi_def_var", ncresult);

  if(! g_mpirank) printf("Writing x...\n");
  ncmpi_enddef(ncid);
  start[0] = g_localx;  count[0] = g_localwidth;
  ncresult = ncmpi_put_vara_double_all(ncid, xvarid, start, count, g_x);
  NCError("ncmpi_put_vara_all (x)", ncresult);
  ncmpi_redef(ncid);

  ncresult = ncmpi_def_dim(ncid, YNAME, g_ysize, &dimid[YDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, YNAME, NC_DOUBLE, 1, &dimid[YDIM], &yvarid);
  NCError("ncmpi_def_var", ncresult);

  if(! g_mpirank) printf("Writing y...\n");
  ncmpi_enddef(ncid);
  start[0] = g_localy;  count[0] = g_localheight;
  ncresult = ncmpi_put_vara_double_all(ncid, yvarid, start, count, g_y);
  NCError("ncmpi_put_vara_all (y)", ncresult);
  ncmpi_redef(ncid);

  ncresult = ncmpi_def_dim(ncid, ZNAME, g_zsize, &dimid[ZDIM]);  
  NCError("ncmpi_def_dim", ncresult);

  ncresult = ncmpi_def_var(ncid, ZNAME, NC_DOUBLE, 1, &dimid[ZDIM], &zvarid);
  NCError("ncmpi_def_var", ncresult);

  if(! g_mpirank) printf("Writing z...\n");
  ncmpi_enddef(ncid);
  ncresult = ncmpi_put_var_double_all(ncid, zvarid, g_z);
  NCError("ncmpi_put_var (z)", ncresult);
  ncmpi_redef(ncid);

  if(! g_mpirank) printf("Defining variables...\n");
  
  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var%02d", i);
    ncresult = ncmpi_def_var(ncid, varname, NC_DOUBLE, NDIMS, dimid, &varid[i]);
    NCError("ncmpi_def_var", ncresult);
  }

  if(! g_mpirank) printf("Closing...\n");

  ncresult = ncmpi_close(ncid);
  NCError("ncmpi_close", ncresult);

  for(i = 0; i < NVARS; i++)
  {
    if(! g_mpirank) printf("Openning...\n");

    ncresult = ncmpi_open(g_mpicomm, filename, NC_WRITE, g_mpiinfo, &ncid);
    NCError("ncmpi_open", ncresult);

    ncmpi_enddef(ncid);
    
    if(! g_mpirank) printf("Writing data...\n");

    start[XDIM] = g_localx; 
    start[YDIM] = g_localy;  
    start[ZDIM] = 0;  
    start[TDIM] = 0;
    
    count[XDIM] = g_localwidth;  
    count[YDIM] = g_localheight;  
    count[ZDIM] = g_zsize;
    count[TDIM] = 1;

    ncresult = ncmpi_put_vara_double_all(ncid, varid[i], start, count, g_data);
    NCError("ncmpi_put_vara_double", ncresult);
    
    if(! g_mpirank) printf("Closing file...\n");
    
    ncresult = ncmpi_close(ncid);
    NCError("ncmpi_close", ncresult);
  }

  elapsed = MPI_Wtime() - starttime;

  if(! g_mpirank) printf("Done.  %f s\n\n", elapsed);

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

void PrintInfo(MPI_Info info, const char* label) 
{
  int nkeys;
  char key[MPI_MAX_INFO_KEY];
  char value[MPI_MAX_INFO_VAL];
  int i;
  int flag;

  if(0 == label || 0 == strlen(label))
    printf("MPI_Info:\n");
  else
    printf("%s:\n", label);

  MPI_Info_get_nkeys(info, &nkeys);

  if(0 == nkeys)
  {
    printf("  None\n");
  }
  else
  {
    for(i = 0; i < nkeys; i++) 
    {
      MPI_Info_get_nthkey(info, i, key);
      MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
      printf("  %s : %s\n", key, value);
    }
  }
}
