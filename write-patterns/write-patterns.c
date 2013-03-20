/**
   write-patterns.c - Created by Timothy Morey on 12/5/2012

   This file defines an MPI program that tests write patterns on a Lustre
   filesystem.  Each process assumes that it has a chunk of data that it must
   write into a shared file.  In the 'contiguous' technique, we will assume that
   the process' data should be contiguous on disk, and in the 'stripe-aligned'
   technique, we will assume that the data consists of noncontiguous stripes of
   contiguous data that happen to be the same size as the Lustre stripes.
 */


#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define USE_MPIIO 1
#define MAX_STRIPES 160

#define MPIERR(code) {                                                  \
    if(MPI_SUCCESS != code) {                                           \
      char errmsg[MPI_MAX_ERROR_STRING];                                \
      int len;                                                          \
      MPI_Error_string(code, errmsg, &len);                             \
      fprintf(stderr, "Rank %03d: MPI ERROR (%d) %s\n",                 \
              mpirank, code, errmsg);                                   \
      exit(code);                                                       \
    }                                                                   \
  }


int main(int argc, char* argv[]) {
  int mpirank, mpisize;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
#ifdef USE_MPIIO
  MPI_File outfile;
  MPI_Status status;
#else
  int outfile;
#endif
  void* buf = 0;
  MPI_Offset buflen = 1073741824; // = 2^30
  int i;
  char outpath[256] = "./outfile";
  char filename[256];
  char value[32];
  int result = 0;
  char c;
  int testcontig = 1;
  int teststriped = 1;
  MPI_Offset fileoffset, bufoffset;
  MPI_Offset stripesize = 1048576;
  int stripecount = 128;
  double start, elapsed, maxelapsed;
  MPI_Offset written;
  int ratio;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &mpirank);
  MPI_Comm_size(comm, &mpisize);

  // Default the number of stripes to the number of processes, and keep halving
  // it until we get a valid stripe count:
  stripecount = mpisize;
  while(stripecount > MAX_STRIPES)
    stripecount /= 2;

  // Parse the command line arguments:
  // Note: if both -c and -r are set, the one that appears later in the command
  //       line arguments will override the former.
  while((c = getopt(argc, argv, "c:m:o:r:s:")) != -1) {
    switch(c) {
    case 'c':  // Stripe count
      stripecount = atoi(optarg);
      break;

    case 'm':  // Which tests are we running?
      if(0 == strcmp("both", optarg) || 
         0 == strcmp("b", optarg)) {
        testcontig = 1;
        teststriped = 1;
      } else if(0 == strcmp("contiguous", optarg) ||
                0 == strcmp("c", optarg)) {
        testcontig = 1;
        teststriped = 0;
      } else if(0 == strcmp("stripe-aligned", optarg) ||
                0 == strcmp("s", optarg)) {
        testcontig = 0;
        teststriped = 1;
      } else {
        if(0 == mpirank)
          fprintf(stderr, "Unrecognized method type (-m value): %s\n", optarg);
      }
      
      break;

    case 'o':  // Output file base name
      strcpy(outpath, optarg);
      break;

    case 'r':  // Process-to-OST ratio
      ratio = atoi(optarg);
      stripecount = mpisize / ratio;
      break;

    case 's':  // Stripe size
      stripesize = atoi(optarg);
      break;

    default:
      if(0 == mpirank)
        fprintf(stderr, "Unrecognized argument %c\n", c);
      break;
    }
  }

  // Check for any problematic parameters:
  if(stripecount > MAX_STRIPES) {
    if(0 == mpirank) {
      fprintf(stderr, "Invalid stripe count: %d.\n", stripecount);
    }

    exit(1);
  }

  // Echo the parameter values we're using to stdout:
  if(0 == mpirank) {
    printf("Using '%s' as the base for all output filenames.\n", outpath);
    printf("Local buffer size is %lld bytes.\n", buflen);
    printf("%s the contiguous write pattern.\n", 
           testcontig ? "Testing" : "Not testing");
    printf("%s the stripe-aligned write pattern.\n",
           teststriped ? "Testing" : "Not testing");
    printf("Stripe size: %lld.\n", stripesize);
    printf("Stripe count: %d.\n", stripecount);
  }

  // Set MPI hints to control the stripe size/count:
  result = MPI_Info_create(&info);  MPIERR(result);  MPIERR(result);
  sprintf(value, "%d", stripecount);
  result = MPI_Info_set(info, "striping_factor", value);  MPIERR(result);
  sprintf(value, "%lld", stripesize);
  result = MPI_Info_set(info, "striping_unit", value);  MPIERR(result);
  result = MPI_Info_set(info, "romio_ds_write", "disable");  MPIERR(result);

  // Initialize our local data:
  if(0 == mpirank) 
    printf("Initializing data...\n");
  buf = malloc(buflen);
  memset(buf, mpirank, buflen);

  if(teststriped) {
    if(0 == mpirank) printf("Performing stripe-aligned write...\n");

    start = MPI_Wtime();

    sprintf(filename, "%s.striped", outpath);
    if(0 == mpirank) printf("Opening %s...\n", filename);

#ifdef USE_MPIIO
    result = MPI_File_open(comm, filename, 
                           MPI_MODE_CREATE | MPI_MODE_WRONLY, 
                           info, &outfile);
    MPIERR(result);
#else
    outfile = open(filename, O_WRONLY | O_CREAT);
    if(-1 == outfile) {
      fprintf(stderr, "Rank %03d: Failed to open '%s' with code %d.\n",
              mpirank, filename, errno);
    }
#endif

    if(0 == mpirank) printf("Writing data...\n");
    // assuming that buflen % stripesize == 0:
    for(i = 0; i < buflen / stripesize; i++) {
      fileoffset = i * mpisize * stripesize + mpirank * stripesize;
      bufoffset = i * stripesize;

#ifdef USE_MPIIO
      result = MPI_File_seek(outfile, fileoffset, MPI_SEEK_SET);
      MPIERR(result);

      result = MPI_File_write(outfile, ((char*)buf) + bufoffset, stripesize, 
                              MPI_BYTE, &status);
      MPIERR(result);
#else
      lseek(outfile, fileoffset, SEEK_SET);
      written = 0;
      while(written < stripesize) {
        result = write(outfile, ((char*)buf) + written, stripesize - written);
        if(result > 0) {
          written += result;
        } else {
          fprintf(stderr, "Rank %03d: write failed with code %d.\n", mpirank, errno);
          break;
        }
      }

      if(written < stripesize) break;
#endif
    }

    if(0 == mpirank) printf("Closing file...\n");

#ifdef USE_MPIIO
    MPI_File_close(&outfile);
#else
    close(outfile);
#endif

    elapsed = MPI_Wtime() - start;

    MPI_Reduce(&elapsed, &maxelapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if(0 == mpirank) printf("Stripe-aligned write took %f s.\n", maxelapsed);
  }

  if(testcontig) {
    if(0 == mpirank) printf("Performing contiguous write...\n");

    start = MPI_Wtime();

    sprintf(filename, "%s.contiguous", outpath);
    if(0 == mpirank) printf("Opening %s...\n", filename);

#ifdef USE_MPIIO
    result = MPI_File_open(comm, filename, 
                           MPI_MODE_CREATE | MPI_MODE_WRONLY, 
                           info, &outfile);
    MPIERR(result);
#else
    outfile = open(filename, O_WRONLY | O_CREAT);
    if(-1 == outfile) {
      fprintf(stderr, "Rank %03d: Failed to open '%s' with code %d.\n",
              mpirank, filename, errno);
    }
#endif

    if(0 == mpirank) printf("Writing data...\n");
    fileoffset = mpirank * buflen;
    
#ifdef USE_MPIIO
    result = MPI_File_seek(outfile, fileoffset, MPI_SEEK_SET);  
    MPIERR(result);
    result = MPI_File_write(outfile, buf, buflen, MPI_BYTE, &status);  
    MPIERR(result);
#else
    lseek(outfile, fileoffset, SEEK_SET);
    written = 0;
    while(written < buflen) {
      result = write(outfile, ((char*)buf) + written, buflen - written);
      if(result > 0) {
        written += result;
      } else {
        fprintf(stderr, "Rand %03d: write failed with code %d.\n", mpirank, errno);
        break;
      }
    }
#endif

    if(0 == mpirank) printf("Closing file...\n");

#ifdef USE_MPIIO
    MPI_File_close(&outfile);
#else
    close(outfile);
#endif

    elapsed = MPI_Wtime() - start;

    MPI_Reduce(&elapsed, &maxelapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if(0 == mpirank) printf("Contiguous write took %f s.\n", maxelapsed);
  }

  MPI_Finalize();

  return 0;
}
