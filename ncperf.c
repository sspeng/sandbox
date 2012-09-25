/**
   ncperf.c - Created by Timothy Morey on 9/22/1012.

   This file defines a program that quantifies the performance of Unidata's 
   NetCDF library with respect to writing out the type of data that PISM tends
   to write.

   Large portions of this file, and certainly inspiration, we borrowed from the
   tst_nc4perf.c file provided with Unidata's NetCDF source package, v4.2.  It 
   contained this descriptive message:

     "This program tests netcdf-4 parallel I/O. These tests are based on the
      needs of the NASA GMAO model, and are based on some test code from
      Dennis Nadeau."
 */

#include <getopt.h>
#include <math.h>
#include <netcdf_par.h> // must be included before netcdf.h
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_MODES 2
#define NUM_FACC 2
#define NUM_CACHE_SIZES 3
#define NUM_CHUNK_COMBOS 3
#define MEGABYTE 1048576
#define XDIM 0
#define YDIM 1
#define TDIM 2
#define MAXVARS 32
#define ROW_FORMAT "%9s | %11s | %10s | %16s | %16s | %16s\n"

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

int TestWrite(size_t cachesize, int facctype, int accessflag, size_t* chunksize,
	      MPI_Comm comm, MPI_Info info, int mpisize, int mpirank, 
	      float* data, int x, int y, int width, int height, int tlen, 
	      int nvars, int ntries, int gwidth, int gheight)
{
  int retval = 0;
  size_t nelems_in;
  float preemption_in;
  int ncid;
  int varid[MAXVARS];
  int ncresult = NC_NOERR;
  int dimids[3];
  size_t start[3];
  size_t count[3];
  double starttime = 0.0;
  double writetime = 0.0;
  double bandwidth = 0.0;
  int try, t, v;
  char filename[NC_MAX_NAME + 1];
  char varname[NC_MAX_NAME + 1];

  // Set the cache size
  if((ncresult = nc_get_chunk_cache(NULL, &nelems_in, &preemption_in)) ||
     (ncresult = nc_set_chunk_cache(cachesize, nelems_in, preemption_in)))
  {
    printf("Failed to set cache size: ncresult=%d, mpirank=%d\n",
	   ncresult, mpirank);
    retval = 1;
  }

  for(try = 0; try < ntries; try++)
  {  
    sprintf(filename, "outfile.t%02d.nc", try);
    // Create a netcdf-4 file, opened for parallel I/O
    if(ncresult = nc_create_par(filename, facctype|NC_NETCDF4, comm, info, &ncid))
    {
      printf("Failed to create file: ncresult=%d, mpirank=%d\n",
	     ncresult, mpirank);
      retval = 1;
    }
    
    // Create our dimensions
    if((ncresult = nc_def_dim(ncid, "x", gwidth, &dimids[XDIM])) ||
       (ncresult = nc_def_dim(ncid, "y", gheight, &dimids[YDIM])) ||
       (ncresult = nc_def_dim(ncid, "t", NC_UNLIMITED, &dimids[TDIM])))
    {
      printf("Failed to create dimensions: ncresult=%d, mpirank=%d\n",
	     ncresult, mpirank);
      retval = 1;
    }
    
    for(v = 0; v < nvars; v++)
    {
      sprintf(varname, "var%02d", v);

      // Create our variables
      if(ncresult = nc_def_var(ncid, varname, NC_FLOAT, 3, dimids, &varid[v]))
      {
	printf("Failed to create variable \"%s\": ncresult=%d, mpirank=%d\n",
	       varname, ncresult, mpirank);
	retval = 1;
      }
      
      // Set the chunk size
      if(chunksize[0])
      {
	if(ncresult = nc_def_var_chunking(ncid, varid[v], 0, chunksize))
	{
	  printf("Failed to set chunk size to %dx%dx%d: ncresult=%d, mpirank=%d\n",
		 chunksize[0], chunksize[1], chunksize[2], ncresult, mpirank);
	  retval = 1;
	}
      }
    }
    
    if(ncresult = nc_enddef(ncid))
    {
      printf("nc_enddef failed: ncresult=%d, mpirank=%d\n",
	     ncresult, mpirank);
      retval = 1;
    }
      
    starttime = MPI_Wtime();

    for(t = 0; t < tlen; t++)
    {    
      start[XDIM] = x;
      start[YDIM] = y;
      start[TDIM] = t;
      
      count[XDIM] = width;
      count[YDIM] = height;
      count[TDIM] = 1;
      
      starttime = MPI_Wtime();

      for(v = 0; v < nvars; v++)
      {
	if(ncresult = nc_var_par_access(ncid, varid[v], accessflag))
	{
	  printf("nc_var_par_access failed: ncresult=%d, mpirank=%d\n",
		 ncresult, mpirank);
	  retval = 1;
	}
	
	// Write as slice of data
	if(ncresult = nc_put_vara_float(ncid, varid[v], start, count, data))
	{
	  printf("nc_put_vara_float failed: ncresult=%d, mpirank=%d\n",
		 ncresult, mpirank);
	  retval = 1;
	}
      }

      // Use a barrier to keep the processes more or less in step, as this will
      // more closely immitate pism, in which all processes are synchronized by
      // time steps.
      MPI_Barrier(MPI_COMM_WORLD);
    }
      
    // Close the file
    if(ncresult = nc_close(ncid))
    {
      printf("nc_close failed: ncresult=%d, mpirank=%d\n",
	     ncresult, mpirank);
      retval = 1;
    }
    
    writetime += MPI_Wtime() - starttime;

    remove(filename);
  }

  writetime /= ntries;
  bandwidth = sizeof(float) * gwidth * gheight * tlen * nvars / MEGABYTE / writetime;
    
  if(0 == mpirank)
  {
    char mpimode_col[32];
    char access_col[32];
    char cache_col[32];
    char chunksize_col[32];
    char writetime_col[32];
    char bandwidth_col[32];
    
    sprintf(mpimode_col, "%s", facctype == NC_MPIIO ? "MPI-IO" : "MPI-POSIX");
    sprintf(access_col, "%s", accessflag == NC_INDEPENDENT ? "independent" : "collective");
    sprintf(cache_col, "%d", (int)cachesize / MEGABYTE);
    if (chunksize[0])
      sprintf(chunksize_col, "%dx%dx%d", 
	      (int)chunksize[0], (int)chunksize[1], (int)chunksize[2]);
    else
      sprintf(chunksize_col, "contiguous");
    sprintf(writetime_col, "%f", writetime);
    sprintf(bandwidth_col, "%f", bandwidth);

    printf(ROW_FORMAT, mpimode_col, access_col, cache_col, chunksize_col, 
	   writetime_col, bandwidth_col);
  }
  
  return retval;
}

int main(int argc, char* argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  int mpisize = 0;
  int mpirank = -1;
  int gwidth = 1501;
  int gheight = 2801;  
  int mpimode[NUM_MODES] = {NC_MPIIO, NC_MPIPOSIX};
  int facctype[NUM_FACC] = {NC_INDEPENDENT, NC_COLLECTIVE};
  size_t cachesize[NUM_CACHE_SIZES] = {MEGABYTE, 32 * MEGABYTE, 64 * MEGABYTE};
  size_t chunksize[NUM_CHUNK_COMBOS][3] = {{0, 0, 0}, 
					   {gwidth, gheight, 1},
					   {32, 32, 1}};
  int m, f, c, i;
  float* data = 0;
  int x, y, width, height;
  int tlen = 10;
  int nvars = 6;
  int ntries = 3;
  char opt;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  while((c = getopt(argc, argv, "h:i:t:v:w:")) != -1)
  {
    switch(c)
    {
    case 'h': // Height of full grid (number of rows)
      gheight = atoi(optarg); break;
    case 'i': // Number of tries that we will average to calculate stats
      ntries = atoi(optarg); break;
    case 't': // Number of time steps
      tlen = atoi(optarg); break;
    case 'v': // Number of variables
      nvars = atoi(optarg); break;
    case 'w': // Width of full grid (number of columns)
      gwidth = atoi(optarg); break;
    }
  }

  chunksize[1][XDIM] = gwidth;
  chunksize[1][YDIM] = gheight;
  chunksize[1][TDIM] = 1;

  // Initialize our local grid and fill it up with simple data
  GetLocalBounds(gwidth, gheight, mpirank, mpisize, &x, &y, &width, &height);
  //printf("Rank %02d: %d %d %d %d\n", mpirank, x, y, width, height);
  data = (float*)malloc(width * height * sizeof(float));
  for(i = 0; i < width*height; i++)
    data[i] = (float)mpirank + ((float)i)/((float)width*height);

  if(0 == mpirank)
  {
    printf("Running tests on %d processes\n", mpisize);
    printf("Write Time and Bandwidth values are averages over %d tries\n", ntries);
    printf("The grid size is %dx%d\n", gwidth, gheight);
    
    for(i = 0; i < mpisize; i++)
    {
      int px, py, pw, ph;
      GetLocalBounds(gwidth, gheight, i, mpisize, &px, &py, &pw, &ph);
      printf("Process %d manages the local grid: x=%4d, y=%4d, width=%4d, height=%4d\n",
	     i, px, py, pw, ph);
    }
    
    printf("Writing %d time slices for %d variables\n", tlen, nvars);
    printf("\n");
    
    printf(ROW_FORMAT, "MPI Mode", "Access", "Cache (MB)", 
	   "Chunk Size", "Write Time (s)", "Bandwidth (MB/s)");
    printf(ROW_FORMAT, "---------", "-----------", "----------", 
	   "----------------", "----------------", "----------------");
  }
  
  for (m = 0; m < NUM_MODES; m++)
    for (f = 0; f < NUM_FACC; f++)
      for (i = 0; i < NUM_CACHE_SIZES; i++)
	for (c = 0; c < NUM_CHUNK_COMBOS; c++)
	  TestWrite(cachesize[i], mpimode[m], facctype[f], chunksize[c], 
		    comm, info, mpisize, mpirank, 
		    data, x, y, width, height, tlen, nvars, ntries, gwidth, gheight);

  return 0;
}
