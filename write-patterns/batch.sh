#!/bin/bash

#$ -V
#$ -cwd
#$ -j y
#$ -o o.$JOB_NAME.$JOB_ID
#$ -q development
#$ -N wp_test
#$ -pe 1way 256
#$ -l h_rt=00:30:00
#$ -M timmorey@gmail.com
#$ -m be

#set -x

OUTDIR=$SCRATCH/$JOB_NAME.$JOB_ID

for i in 1
do
	mkdir $OUTDIR
	ibrun tacc_affinity $WORK/pism/sandbox/write-patterns/write-patterns-mpiio -o $OUTDIR/outfile -m stripe-aligned
	ibrun tacc_affinity $WORK/pism/sandbox/write-patterns/write-patterns-mpiio -o $OUTDIR/outfile -m contiguous
	#rm -rf $OUTDIR
done

