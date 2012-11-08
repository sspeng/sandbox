#!/bin/bash

# Copyright (C) 2011-2012 Andy Aschwanden and Gudfinna Adalgeirsdottir
#
# before using this script, run preprocess.sh to download and adjust metadata
# on input files, then run prespinup.

# seconds per year, from UDUNITS
SECPERA=3.15569259747e7

RUN=1
POST=0
if [ -n "${PISM_RUN:+1}" ] ; then
    RUN=$PISM_RUN
fi
if [ -n "${PISM_POST:+1}" ] ; then
    POST=$PISM_POST
fi


EXPERIMENT=FOO
if [ -n "${PISM_EXPERIMENT:+1}" ] ; then
    EXPERIMENT=$PISM_EXPERIMENT
else
    echo
    echo "EXPERIMENT=$EXPERIMENT"
fi
SCRIPTNAME="#(${EXPERIMENT})"

echo
echo "# =================================================================================="
echo "# PISM initialization run (aka spinup) '$EXPERIMENT'"
echo "# =================================================================================="
echo

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME                              NN = $NN"

# set output format:
#  $ export PISM_OFORMAT="netcdf4_parallel "
if [ -n "${PISM_OFORMAT:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                      PISM_OFORMAT = $PISM_OFORMAT  (already set)"
else
  PISM_OFORMAT="netcdf3"
  echo "$SCRIPTNAME                      PISM_OFORMAT = $PISM_OFORMAT"
fi
OFORMAT=$PISM_OFORMAT

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME                      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME                         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                    PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME                     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME                       PISM_EXEC = $PISM_EXEC"
fi

OCEAN="-ocean_kill $PISM_FTT_FILE -shelf_basal_melt_rate 2"
PISM_INITIALIZATION_BCFILE=ERAI_1989_2011_5KM_CON_MEAN.nc
PISM_FORCING_BCFILE=ERAI_1989_2011_5KM_CON_MM.nc

# set COUPLERS
if [ -n "${PISM_INITIALIZATION_COUPLER:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_INITIALIZATION_COUPLER = $PISM_INITIALIZATION_COUPLER  (already set)"
else
  PISM_INITIALIZATION_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_given_file $PISM_INITIALIZATION_BCFILE -surface_lapse_rate_file $PISM_INITIALIZATION_BCFILE $OCEAN"
  echo "$SCRIPTNAME     PISM_INITIALIZATION_COUPLER = $PISM_INITIALIZATION_COUPLER"
fi
INITIALIZATION_COUPLER=$PISM_INITIALIZATION_COUPLER
if [ -n "${PISM_FTT_COUPLER:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME               PISM_FTT_COUPLER = $PISM_INITIALIZATION_COUPLER  (already set)"
else
  PISM_FTT_COUPLER=$PISM_INITIALIZATION_COUPLER
  echo "$SCRIPTNAME               PISM_FTT_COUPLER = $PISM_FTT_COUPLER"
fi
FTT_COUPLER=$PISM_FTT_COUPLER
if [ -n "${PISM_FORCING_COUPLER:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME            PISM_FORCING_COUPLER = $PISM_FORCING_COUPLER  (already set)"
else
  PISM_FORCING_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_given_file $PISM_FORCING_BCFILE -surface_lapse_rate_file $PISM_FORCING_BCFILE $OCEAN"
  echo "$SCRIPTNAME            PISM_FORCING_COUPLER = $PISM_FORCING_COUPLER"
fi
FORCING_COUPLER=$PISM_FORCING_COUPLER
if [ -n "${PISM_REFERENCE_COUPLER:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME            PISM_REFERENCE_COUPLER = $PISM_REFERENCE_COUPLER  (already set)"
else
  PISM_REFERENCE_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1  -surface_lapse_rate_file $PISM_REFERENCE_BCFILE $OCEAN"
  echo "$SCRIPTNAME            PISM_REFERENCE_COUPLER = $PISM_REFERENCE_COUPLER"
fi
REFERENCE_COUPLER=$PISM_REFERENCE_COUPLER

if [ -n "${PISM_FORCING_BCFILE:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME            PISM_FORCING_BCFILE = $PISM_FORCING_BCFILE  (already set)"
else
  PISM_FORCING_BCFILE="ERAI_1989_2011_5KM_CON_MM.nc"
  echo "$SCRIPTNAME            PISM_FORCING_BCFILE = $PISM_FORCING_BCFILE"
fi

if [ -n "${PISM_TITLE:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                      PISM_TITLE = $PISM_TITLE  (already set)"
else
  PISM_TITLE="Constant-climate initialization with 1989-2011 CMB"
  echo "$SCRIPTNAME                      PISM_TITLE = $PISM_TITLE"
fi
TITLE=$PISM_TITLE

# run length starting time
if [ -n "${PISM_PALEOSTARTYEAR:+1}" ] ; then  # check if env var is already set
  PALEOSTARTYEAR=$PISM_PALEOSTARTYEAR
  echo "$SCRIPTNAME            PISM_PALEOSTARTYEAR = $PISM_PALEOSTARTYEAR  (already set)"
else
  PALEOSTARTYEAR=-125000
  echo "$SCRIPTNAME            PISM_PALEOSTARTYEAR = $PALEOSTARTYEAR"
fi

COARSESTENDTIME=-7000
COARSEENDTIME=-2000
MEDIUMENDTIME=-850
FINEENDTIME=-350
SUPERENDTIME=-100

# grids
GRID20KM="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 21 -z_spacing equal"
GRID10KM="-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 201 -Mbz 21 -z_spacing equal"
GRID5KM="-Mx 301 -My 561 -Lz 4000 -Lbz 2000 -Mz 201 -Mbz 21 -z_spacing equal"
GRID3KM="-Mx 501 -My 935 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID2FKM="-Mx 601 -My 1121 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID2KM="-Mx 751 -My 1401 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID1KM="-Mx 1501 -My 2801 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"


# grid spacings
GS20KM=20
GS10KM=10
GS5KM=5
GS3KM=3
GS2FKM=2.5
GS2KM=2
GS1KM=1

# skips
SKIP20KM=10
SKIP10KM=50
SKIP5KM=200
SKIP3KM=500
SKIP2FKM=750
SKIP2KM=1000
SKIP1KM=2000

# defaults to coarse grid choices
COARSESTGRID=$GRID20KM
COARSEGRID=$GRID20KM
MEDIUMGRID=$GRID20KM
FINEGRID=$GRID20KM
SUPERGRID=$GRID20KM
ULTRAGRID=$GRID20KM
COARSESTSKIP=$SKIP20KM
COARSESKIP=$SKIP20KM
MEDIUMSKIP=$SKIP20KM
FINESKIP=$SKIP20KM
SUPERSKIP=$SKIP20KM
ULTRASKIP=$SKIP20KM
VCGS=$GS20KM
CGS=$GS20KM
MGS=$GS20KM
FGS=$GS20KM
SGS=$GS20KM
UGS=$GS20KM

echo ""
if [ $# -gt 1 ] ; then
  if [ $2 -eq "1" ] ; then  # if user says "spinup.sh N 1" then use:
    echo "$SCRIPTNAME grid: RUNS ON 20km and 10km"
    echo "$SCRIPTNAME       WARNING: MEDIUM COMPUTATIONAL TIME"
    COARSESTGRID=$GRID20KM
    COARSEGRID=$GRID10KM
    MEDIUMGRID=$GRID10KM
    FINEGRID=$GRID10KM
    SUPERGRID=$GRID10KM
    ULTRAGRID=$GRID10KM
    COARSESTSKIP=$SKIP20KM
    COARSESKIP=$SKIP10KM
    MEDIUMSKIP=$SKIP10KM
    FINESKIP=$SKIP10KM
    SUPERSKIP=$SKIP10KM
    ULTRASKIP=$SKIP10KM
    VCGS=$GS20KM
    CGS=$GS10KM
    MGS=$GS10KM
    FGS=$GS10KM
    SGS=$GS10KM
    UGS=$GS10KM
  fi
  if [ $2 -eq "2" ] ; then  # if user says "spinup.sh N 2" then use:
    echo "$SCRIPTNAME grid: RUNS ON 20km, 10km, 5km and 5km"
    echo "$SCRIPTNAME       WARNING: LARGE COMPUTATIONAL TIME"
    COARSESTGRID=$GRID20KM
    COARSEGRID=$GRID10KM
    MEDIUMGRID=$GRID5KM
    FINEGRID=$GRID5KM
    SUPERGRID=$GRID5KM
    ULTRAGRID=$GRID5KM
    COARSESTSKIP=$SKIP20KM
    COARSESKIP=$SKIP10KM
    MEDIUMSKIP=$SKIP5KM
    FINESKIP=$SKIP5KM
    SUPERSKIP=$SKIP5KM
    ULTRASKIP=$SKIP5KM
    VCGS=$GS20KM
    CGS=$GS10KM
    MGS=$GS5KM
    FGS=$GS5KM
    SGS=$GS5KM
    UGS=$GS5KM
  fi
  if [ $2 -eq "3" ] ; then  # if user says "spinup.sh N 3":
    echo "$SCRIPTNAME grid: RUNS ON 20km, 10km, 5km and 2.5km"
    echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
    COARSESTGRID=$GRID20KM
    COARSEGRID=$GRID10KM
    MEDIUMGRID=$GRID5KM
    FINEGRID=$GRID2FKM
    SUPERGRID=$GRID2FKM
    ULTRAGRID=$GRID2FKM
    COARSESTSKIP=$SKIP20KM
    COARSESKIP=$SKIP10KM
    MEDIUMSKIP=$SKIP5KM
    FINESKIP=$SKIP2FKM
    SUPERSKIP=$SKIP2FKM
    ULTRASKIP=$SKIP2FKM
    VCGS=$GS20KM
    CGS=$GS10KM
    MGS=$GS5KM
    FGS=$GS2FKM
    SGS=$GS2FKM
    UGS=$GS2FKM
  fi
  if [ $2 -eq "4" ] ; then  # if user says "spinup.sh N 4":
    echo "$SCRIPTNAME grid: RUNS ON 10km, 5km, 2.5km and 2km"
    echo "$SCRIPTNAME       WARNING: EXTREMELY LARGE COMPUTATIONAL TIME"
    COARSESTGRID=$GRID20KM
    COARSEGRID=$GRID10KM
    MEDIUMGRID=$GRID5KM
    FINEGRID=$GRID2FKM
    SUPERGRID=$GRID2KM
    ULTRAGRID=$GRID2KM
    COARSESTSKIP=$SKIP20KM
    COARSESKIP=$SKIP10KM
    MEDIUMSKIP=$SKIP5KM
    FINESKIP=$SKIP2FKM
    SUPERSKIP=$SKIP2KM
    ULTRASKIP=$SKIP2KM
    VCGS=$GS20KM
    CGS=$GS10KM
    MGS=$GS5KM
    FGS=$GS2FKM
    SGS=$GS2KM
    UGS=$GS2KM
  fi
  if [ $2 -eq "5" ] ; then  # if user says "spinup.sh N 5":
    echo "$SCRIPTNAME grid: RUNS ON 10km, 5km, 2.5km, 2km and 1km"
    echo "$SCRIPTNAME       WARNING: EXTREMELY LARGE COMPUTATIONAL TIME"
    COARSESTGRID=$GRID20KM
    COARSEGRID=$GRID10KM
    MEDIUMGRID=$GRID5KM
    FINEGRID=$GRID2FKM
    SUPERGRID=$GRID2KM
    ULTRAGRID=$GRID1KM
    COARSESTSKIP=$SKIP20KM
    COARSESKIP=$SKIP10KM
    MEDIUMSKIP=$SKIP5KM
    FINESKIP=$SKIP2FKM
    SUPERSKIP=$SKIP2KM
    ULTRASKIP=$SKIP1KM
    VCGS=$GS20KM
    CGS=$GS10KM
    MGS=$GS5KM
    FGS=$GS2FKM
    SGS=$GS2KM
    UGS=$GS1KM
  fi
else
    echo "$SCRIPTNAME grid: ALL RUNS ON 20km"
fi
echo ""

# preprocess.sh generates pism_*.nc files; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATANAME=pism_Greenland_${UGS}km_griggs_updated_smoothed.nc
fi

INNAME=$PISM_DATANAME


INITIAL_CONFIG=initial_config.nc
ERA_CONFIG=era_config.nc
STEADYNAME=g${CGS}km_steady_ssa.nc
for INPUT in $PISM_DATANAME $CONFIG $STEADYNAME; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "     !!!!   RUN  preprocess.sh  TO GENERATE  $INPUT   !!!!"
    echo
  fi
done

echo "$SCRIPTNAME   coarsest grid = '$COARSESTGRID' (= $VCGS km)"
echo "$SCRIPTNAME     coarse grid = '$COARSEGRID' (= $CGS km)"
echo "$SCRIPTNAME     medium grid = '$MEDIUMGRID' (= $MGS km)"
echo "$SCRIPTNAME       fine grid = '$FINEGRID' (= $FGS km)"
echo "$SCRIPTNAME      finer grid = '$SUPERGRID' (= $SGS km)"
echo "$SCRIPTNAME     finest grid = '$ULTRAGRID' (= $UGS km)"
# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -prof -title \"$TITLE\" "

# default choices in parameter study; see Bueler & Brown (2009) re "tillphi"
TILLPHI="-topg_to_phi 5.0,20.0,-300.0,700.0"
FULLPHYS="-ssa_sliding -thk_eff $TILLPHI"

echo "$SCRIPTNAME             executable = '$PISM'"
echo "$SCRIPTNAME           full physics = '$FULLPHYS'"
echo "$SCRIPTNAME initialization coupler = '$INITIALIZATION_COUPLER'"
echo "$SCRIPTNAME            ftt coupler = '$FTT_COUPLER'"
echo "$SCRIPTNAME        forcing coupler = '$FORCING_COUPLER'"
echo "$SCRIPTNAME         future coupler = '$FUTURE_COUPLER'"


# time step for scalars
TSSTEP=yearly
# time step for exvars to save
EXSTEP=500
# extra variables saved every EXSTEP
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask,cbase,tauc,thk,Href,IcebergMask,mask,lat,lon,climatic_mass_balance_cumulative,ocean_kill_flux_cumulative,taud_mag,topg,usurf,velbar"
# regrid vars
REGRIDVARS="litho_temp,thk,enthalpy,bwat,bmelt"
# output file size
OSIZE="medium"

# start back at 125ka BPE
STARTTIME=$PALEOSTARTYEAR
ENDTIME=$COARSESTENDTIME # BP
ET=$(($ENDTIME/-1))
STARTNAME=$STEADYNAME
OUTNAME=g${VCGS}km_m${ET}a_$EXPERIMENT.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME

cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $COARSESKIP -boot_file $INNAME $COARSESTGRID $FULLPHYS \
     $INITIALIZATION_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  run with full physics,"
    echo "$SCRIPTNAME    from ${PALEOSTARTYEAR}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi

STARTTIME=$ENDTIME
ENDTIME=$COARSEENDTIME # BP
ET=$(($ENDTIME/-1))
STARTNAME=$OUTNAME
OUTNAME=g${CGS}km_m${ET}a_$EXPERIMENT.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=100
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME

cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $MEDIUMSKIP -boot_file $INNAME $COARSEGRID $FULLPHYS \
     $INITIALIZATION_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  regrid and do run with full physics,"
    echo "$SCRIPTNAME    from ${STARTTIME}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi

STARTTIME=$ENDTIME
ENDTIME=$MEDIUMENDTIME # BP
ET=$(($ENDTIME/-1))
STARTNAME=$OUTNAME
OUTNAME=g${MGS}km_m${ET}a_$EXPERIMENT.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=100
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME

cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $MEDIUMSKIP -boot_file $INNAME $MEDIUMGRID $FULLPHYS \
      $FTT_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  regrid and do run with full physics,"
    echo "$SCRIPTNAME    from ${STARTTIME}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi


STARTTIME=$ENDTIME
ENDTIME=$FINEENDTIME # BP
ET=$(($ENDTIME/-1))
STARTNAME=$OUTNAME
OUTNAME=g${FGS}km_m${ET}a_$EXPERIMENT.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=10
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $FINESKIP -boot_file $INNAME $FINEGRID $FULLPHYS \
      $FTT_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  regrid and do run with full physics,"
    echo "$SCRIPTNAME    from ${STARTTIME}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi


STARTTIME=$ENDTIME
ENDTIME=$SUPERENDTIME # BP
ET=$(($ENDTIME/-1))
STARTNAME=$OUTNAME
OUTNAME=g${SGS}km_m${ET}a_$EXPERIMENT.nc
TSNAME=ts_$OUTNAME
# time step for scalars
TSSTEP=yearly
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=20
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME

cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $SUPERSKIP -boot_file $INNAME $SUPERGRID $FULLPHYS \
      $FTT_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  regrid and do run with full physics,"
    echo "$SCRIPTNAME    from ${STARTTIME}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi


STARTTIME=$ENDTIME
ENDTIME=0 # BP
STARTNAME=$OUTNAME
OUTNAME=g${UGS}km_0_$EXPERIMENT.nc
START1KM=$OUTNAME
TSNAME=ts_$OUTNAME
# time step for scalars
TSSTEP=yearly
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=10
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
LASTNAME=$OUTNAME
cmd="$PISM_MPIDO $NN $PISM -config_override $INITIAL_CONFIG -skip -skip_max $ULTRASKIP -boot_file $INNAME $ULTRAGRID $FULLPHYS \
      $FTT_COUPLER \
     -regrid_file $STARTNAME -regrid_vars $REGRIDVARS   \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
if [ $RUN -eq 1 ]; then
    echo
    echo "$SCRIPTNAME  regrid and do run with full physics,"
    echo "$SCRIPTNAME    from ${STARTTIME}a BPE to ${ENDTIME}a BPE"

    $PISM_DO $cmd
fi

echo


if [ $POST -eq 1 ]; then
    echo "$SCRIPTNAME  some postprocessing"
     # calculate averages of climatic_mass_balance and dHdt using ncap2 sleight of hand.
    cmd="ncap2 -O -s '*sz_idt=time.size(); *climatic_mass_balance[\$time,\$x,\$y]= 0.f; *dHdt[\$time,\$x,\$y]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {climatic_mass_balance(idt,:,:)=(climatic_mass_balance_cumulative(idt,:,:)-climatic_mass_balance_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA; dHdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA;} climatic_mass_balance.ram_write(); dHdt.ram_write();' $EXNAME $EXNAME"
    $PISM_DO $cmd
    echo
     # adjust meta data for new fields
    cmd="ncatted -a units,climatic_mass_balance,o,c,'m year-1' -a units,dHdt,o,c,'m year-1' \
      -a long_name,climatic_mass_balance,o,c,'surface mass balance' \
      -a long_name,dHdt,o,c,'rate of change of ice thickness' \
      -a grid_mapping,climatic_mass_balance,o,c,'mapping' \
      -a grid_mapping,dHdt,o,c,'mapping' \
      -a cell_methods,climatic_mass_balance,o,c,'time: mean (interval: 10 years)' \
      -a cell_methods,dHdt,o,c,'time: mean (interval: 10 years)' $EXNAME"
    $PISM_DO $cmd
    echo
    # adjust meta data for new fields
    cmd="ncks -A -v climatic_mass_balance -d time,0. $EXNAME $OUTNAME"
    $PISM_DO $cmd
    cmd="ncpdq -O -4 -L 3 -a time,y,x $EXNAME ../transfer/$EXNAME"
    echo
    $PISM_DO $cmd
    cmd="ncpdq -O -4 -L 3 -a time,z,y,x $OUTNAME ../transfer/$OUTNAME"
    $PISM_DO $cmd

fi

echo
echo "$SCRIPTNAME  initialization done"


STARTNAME=$OUTNAME
RUNLENGTH=23
TSSTEP=daily
EXSTEP=monthly    
EXVARS="ocean_kill_flux_cumulative,climatic_mass_balance_cumulative,ocean_kill_flux_cumulative,diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask,cbase,tauc,thk,mask,usurf,lat,lon,velbar"
# for STARTTIME in 0 2 4 6; do
for STARTTIME in 0; do

    STARTYEAR=$(($STARTTIME + 1989))
    ENDYEAR=$(($STARTYEAR + $RUNLENGTH - 1))

    echo
    echo "# ============================================================================="
    echo "# PISM HIRHAM Greenland: $STARTYEAR-$ENDYEAR forcing runs"
    echo "# ============================================================================="
    echo
    
    OUTNAME=g${UGS}km_${STARTYEAR}-${ENDYEAR}_$EXPERIMENT.nc
    TSNAME=ts_$OUTNAME
    TSTIMES=$TSSTEP
    EXNAME=ex_$OUTNAME
    EXTIMES=$EXSTEP
    echo
    echo "$SCRIPTNAME  forcing run"
    cmd="$PISM_MPIDO $NN $PISM -config_override $ERA_CONFIG -climatic_mass_balance_cumulative -skip -skip_max $ULTRASKIP -i $STARTNAME $FULLPHYS \
      $FORCING_COUPLER -calendar gregorian -time_file $PISM_FORCING_BCFILE \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
    if [ $RUN -eq 1 ]; then
        $PISM_DO $cmd
    fi

    if [ $POST -eq 1 ]; then
        echo
        echo "# ============================================================================="
        echo "# PISM HIRHAM Greenland: $STARTYEAR-$ENDYEAR forcing runs"
        echo "# ============================================================================="
        echo

        # calculate climatic_mass_balance
        echo
        echo "$SCRIPTNAME  postprocessing forcing run"
        # calculate monthly-averages of climatic_mass_balance and dHdt using ncap2 sleight of hand.
	cmd="ncap2 -O -s '*sz_idt=time.size(); *climatic_mass_balance[\$time,\$x,\$y]= 0.f; *dHdt[\$time,\$x,\$y]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {climatic_mass_balance(idt,:,:)=(climatic_mass_balance_cumulative(idt,:,:)-climatic_mass_balance_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA; dHdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA;} climatic_mass_balance.ram_write(); dHdt.ram_write();' $EXNAME $EXNAME"
        $PISM_DO $cmd
	echo
       # adjust meta data for new fields
       cmd="ncatted -a units,climatic_mass_balance,o,c,'m year-1' \
        -a units,dHdt,o,c,'m year-1' \
        -a long_name,climatic_mass_balance,o,c,'surface mass balance' \
        -a long_name,dHdt,o,c,'rate of change of ice thickness' \
        -a grid_mapping,climatic_mass_balance,o,c,'mapping' \
        -a grid_mapping,dHdt,o,c,'mapping' \
        -a cell_methods,climatic_mass_balance,o,c,'time: mean (interval: 1 month)' \
        -a cell_methods,dHdt,o,c,'time: mean (interval: 1 month)' $EXNAME"

        $PISM_DO $cmd
	cmd="ncpdq -O -4 -L 3 -a time,y,x $EXNAME ../transfer/$EXNAME"
	echo
        $PISM_DO $cmd
	cmd="ncpdq -O -4 -L 3 -a time,z,y,x $OUTNAME ../transfer/$OUTNAME"
        $PISM_DO $cmd
	cmd="ncks -O $TSNAME ../transfer/$TSNAME"
        $PISM_DO $cmd
    fi

    echo
    echo "# ============================================================================="
    echo "# PISM HIRHAM Greenland: $STARTYEAR-$ENDYEAR reference runs"
    echo "# ============================================================================="
    echo
    
    OUTNAME=g${UGS}km_${STARTYEAR}-${ENDYEAR}_${EXPERIMENT}_REF.nc
    TSNAME=ts_$OUTNAME
    TSTIMES=$TSSTEP
    EXNAME=ex_$OUTNAME
    EXTIMES=$EXSTEP
    echo
    echo "$SCRIPTNAME  forcing run"
    cmd="$PISM_MPIDO $NN $PISM -config_override $ERA_CONFIG -skip -skip_max $ULTRASKIP -i $STARTNAME $FULLPHYS \
      $REFERENCE_COUPLER -calendar gregorian -time_file $PISM_FORCING_BCFILE \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT 2>&1 | tee job.\${PBS_JOBID}"
    if [ $RUN -eq 1 ]; then
        $PISM_DO $cmd
    fi

    if [ $POST -eq 1 ]; then
        echo
        echo "# ============================================================================="
        echo "# PISM HIRHAM Greenland: $STARTYEAR-$ENDYEAR reference runs"
        echo "# ============================================================================="
        echo

        # calculate climatic_mass_balance
        echo
        echo "$SCRIPTNAME  postprocessing forcing run"
        # calculate monthly-averages of climatic_mass_balance and dHdt using ncap2 sleight of hand.
	cmd="ncap2 -O -s '*sz_idt=time.size(); *climatic_mass_balance[\$time,\$x,\$y]= 0.f; *dHdt[\$time,\$x,\$y]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {climatic_mass_balance(idt,:,:)=(climatic_mass_balance_cumulative(idt,:,:)-climatic_mass_balance_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA; dHdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA;} climatic_mass_balance.ram_write(); dHdt.ram_write();' $EXNAME $EXNAME"
        $PISM_DO $cmd
	echo
       # adjust meta data for new fields
       cmd="ncatted -a units,climatic_mass_balance,o,c,'m year-1' \
        -a units,dHdt,o,c,'m year-1' \
        -a long_name,climatic_mass_balance,o,c,'surface mass balance' \
        -a long_name,dHdt,o,c,'rate of change of ice thickness' \
        -a grid_mapping,climatic_mass_balance,o,c,'mapping' \
        -a grid_mapping,dHdt,o,c,'mapping' \
        -a cell_methods,climatic_mass_balance,o,c,'time: mean (interval: 1 month)' \
        -a cell_methods,dHdt,o,c,'time: mean (interval: 1 month)' $EXNAME"

        $PISM_DO $cmd
	cmd="ncpdq -O -4 -L 3 -a time,y,x $EXNAME ../transfer/$EXNAME"
	echo
        $PISM_DO $cmd
	cmd="ncpdq -O -4 -L 3 -a time,z,y,x $OUTNAME ../transfer/$OUTNAME"
        $PISM_DO $cmd
	cmd="ncks -O $TSNAME ../transfer/$TSNAME"
        $PISM_DO $cmd

    fi
    echo
done

echo
echo "$SCRIPTNAME  ERA-interim done"