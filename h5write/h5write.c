/**
   h5write.c - Created by Timothy Morey on 2/22/2013

   This file attempts to replicate PISM's write pattern to test the performance
   of the hdf5 output method.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

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

#define OLD_WRITE_PATTERN 1


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
  const char* filename = "output.h5";
  hid_t fileid, plist, filespace, memspace, dimvar, varid;
  hsize_t size[NDIMS], maxsize[NDIMS], chunksize[NDIMS];
  hsize_t start[NDIMS], count[NDIMS];
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

  if(! mpirank) printf("Creating HDF5 file...\n");

  plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, mpicomm, mpiinfo);

  // TODO: this seems like a good place to put optimizations, and indeed
  // PISM is adding several additional properties, like setting block sizes,
  // cache eviction policies, fs striping parameters, etc.

  fileid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  H5Pclose(plist);

  if(! mpirank) printf("Setting up dimensions...\n");

  if(! mpirank) printf("Creating time dimension...\n");

  // Define the time dimension
  size[0] = 1;
  maxsize[0] = H5S_UNLIMITED;
  chunksize[0] = 1;

  filespace = H5Screate_simple(1, size, maxsize);
  plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist, 1, chunksize);  // It is strictly required to set chunksize when using
                                      // the low-level api.  Contiguous datasets are not allowed
                                      // to use the unlimited dimension.
  dimvar = H5Dcreate(fileid, TNAME, H5T_NATIVE_DOUBLE, filespace, 
                     H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Pclose(plist);
  H5DSset_scale(dimvar, TNAME);
  H5Dclose(dimvar);
  H5Sclose(filespace);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing time dimension...\n");
  dimvar = H5Dopen(fileid, TNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, size, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE); // TODO: Pism does this, but comments suggest it is questionable
  start[0] = 0;
  count[0] = 1;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, &t);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);
#endif

  if(! mpirank) printf("Creating x dimension...\n");

  size[0] = XSIZE;
  
  filespace = H5Screate_simple(1, size, 0);
  dimvar = H5Dcreate(fileid, XNAME, H5T_NATIVE_DOUBLE, filespace,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5DSset_scale(dimvar, XNAME);
  H5Dclose(dimvar);
  H5Sclose(filespace);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing x dimension...\n");
  dimvar = H5Dopen(fileid, XNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, size, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  start[0] = 0;
  count[0] = XSIZE;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, x);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);
#endif

  if(! mpirank) printf("Creating y dimension...\n");

  size[0] = YSIZE;
  
  filespace = H5Screate_simple(1, size, 0);
  dimvar = H5Dcreate(fileid, YNAME, H5T_NATIVE_DOUBLE, filespace,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5DSset_scale(dimvar, YNAME);
  H5Dclose(dimvar);
  H5Sclose(filespace);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing y dimension...\n");
  dimvar = H5Dopen(fileid, YNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, size, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  start[0] = 0;
  count[0] = YSIZE;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, y);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);
#endif

  if(! mpirank) printf("Creating z dimension...\n");

  size[0] = ZSIZE;
  
  filespace = H5Screate_simple(1, size, 0);
  dimvar = H5Dcreate(fileid, ZNAME, H5T_NATIVE_DOUBLE, filespace,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5DSset_scale(dimvar, ZNAME);
  H5Dclose(dimvar);
  H5Sclose(filespace);

#ifdef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing z dimension...\n");
  dimvar = H5Dopen(fileid, ZNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, size, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  start[0] = 0;
  count[0] = ZSIZE;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, z);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);
#endif

  if(! mpirank) printf("Defining variables...\n");

  MPI_Allreduce(&localwidth, &maxwidth, 1, MPI_INT, MPI_MAX, mpicomm);
  MPI_Allreduce(&localheight, &maxheight, 1, MPI_INT, MPI_MAX, mpicomm);
 
  size[TDIM] = 1;
  size[XDIM] = XSIZE;
  size[YDIM] = YSIZE;
  size[ZDIM] = ZSIZE;

  maxsize[TDIM] = H5S_UNLIMITED;
  maxsize[XDIM] = XSIZE;
  maxsize[YDIM] = YSIZE;
  maxsize[ZDIM] = ZSIZE;

  chunksize[TDIM] = 1;
  chunksize[XDIM] = maxwidth;
  chunksize[YDIM] = maxheight;
  chunksize[ZDIM] = ZSIZE;  // Looks like pism might use 1 here...

  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var3d-%02d", i);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, NDIMS, chunksize);
    filespace = H5Screate_simple(NDIMS, size, maxsize);
    varid = H5Dcreate(fileid, varname, H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(filespace);
    H5Dclose(varid);

    sprintf(varname, "var2d-%02d", i);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, NDIMS-1, chunksize);
    filespace = H5Screate_simple(NDIMS-1, size, maxsize);
    varid = H5Dcreate(fileid, varname, H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(filespace);
    H5Dclose(varid);
  }

#ifndef OLD_WRITE_PATTERN
  if(! mpirank) printf("Writing time dimension...\n");
  start[0] = 0;
  count[0] = 1;
  dimvar = H5Dopen(fileid, TNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, count, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE); // TODO: Pism does this, but comments suggest it is questionable
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, &t);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);

  if(! mpirank) printf("Writing x dimension...\n");
  start[0] = 0;
  count[0] = XSIZE;
  dimvar = H5Dopen(fileid, XNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, count, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, x);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);

  if(! mpirank) printf("Writing y dimension...\n");
  start[0] = 0;
  count[0] = YSIZE;
  dimvar = H5Dopen(fileid, YNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, count, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, y);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);

  if(! mpirank) printf("Writing z dimension...\n");
  start[0] = 0;
  count[0] = ZSIZE;
  dimvar = H5Dopen(fileid, ZNAME, H5P_DEFAULT);
  filespace = H5Dget_space(dimvar);
  memspace = H5Screate_simple(1, count, 0);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
  H5Dwrite(dimvar, H5T_NATIVE_DOUBLE, memspace, filespace, plist, z);
  H5Pclose(plist);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Dclose(dimvar);
#endif

  if(! mpirank) printf("Writing variable data...\n");

  for(i = 0; i < NVARS; i++)
  {
    sprintf(varname, "var3d-%02d", i);
    if(! mpirank) printf("Writing %s...\n", varname);
    start[TDIM] = 0;
    start[XDIM] = localx;
    start[YDIM] = localy;
    start[ZDIM] = 0;
    count[TDIM] = 1;
    count[XDIM] = localwidth;
    count[YDIM] = localheight;
    count[ZDIM] = ZSIZE;
    varid = H5Dopen(fileid, varname, H5P_DEFAULT);
    filespace = H5Dget_space(varid);
    memspace = H5Screate_simple(NDIMS, count, 0);
    plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
    H5Dwrite(varid, H5T_NATIVE_DOUBLE, memspace, filespace, plist, data3d);
    H5Pclose(plist);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Dclose(varid);

    sprintf(varname, "var2d-%02d", i);
    if(! mpirank) printf("Writing %s...\n", varname);
    start[TDIM] = 0;
    start[XDIM] = localx;
    start[YDIM] = localy;
    count[TDIM] = 1;
    count[XDIM] = localwidth;
    count[YDIM] = localheight;
    varid = H5Dopen(fileid, varname, H5P_DEFAULT);
    filespace = H5Dget_space(varid);
    memspace = H5Screate_simple(NDIMS-1, count, 0);
    plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0, count, 0);
    H5Dwrite(varid, H5T_NATIVE_DOUBLE, memspace, filespace, plist, data2d);
    H5Pclose(plist);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Dclose(varid);
  }

  if(! mpirank) printf("Closing file...\n");

  H5Fclose(fileid);

  if(! mpirank) printf("Done.\n");

  free(data2d);
  free(data3d);
  free(x);
  free(y);
  free(z);

  MPI_Finalize();

  return 0;
}
