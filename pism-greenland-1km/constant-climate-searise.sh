#!/bin/bash

# Copyright (C) 2011-2012 Andy Aschwanden
#
# *****************************************************************************
# Constant-climate initialization
# *****************************************************************************
#
# - using 1989-2009 monthly surface mass balance

PISM_EXPERIMENT=CONST
PISM_TITLE="Constant-climate initialization"
PISM_INITIALIZATION_BCFILE=pism_Greenland_1km_griggs_updated_smoothed.nc
PISM_FORCING_BCFILE=pism_Greenland_1km_griggs_updated_smoothed.nc
PISM_OK_FILE=pism_Greenland_1km_griggs_updated_smoothed.nc

for INPUT in $PISM_INITIALIZATION_BCFILE $PISM_FORCING_BCFILE $PISM_FTT_FILE; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "# input   $INPUT (found)"
  else
    echo "# input   $INPUT (MISSING!!)"
    exit
  fi
done

OCEAN="-ocean_kill $PISM_OK_FILE -shelf_basal_melt_rate 2"
PISM_INITIALIZATION_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_given_file $PISM_INITIALIZATION_BCFILE -surface_lapse_rate_file $PISM_INITIALIZATION_BCFILE $OCEAN"
PISM_FTT_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_given_file $PISM_INITIALIZATION_BCFILE -surface_lapse_rate_file $PISM_INITIALIZATION_BCFILE $OCEAN"
PISM_FORCING_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_given_file $PISM_FORCING_BCFILE -surface_lapse_rate_file $PISM_FORCING_BCFILE $OCEAN"
PISM_REFERENCE_COUPLER="-surface given,lapse_rate -temp_lapse_rate 7.1 -surface_lapse_rate_file $PISM_FORCING_BCFILE $OCEAN"
PISM_PALEOSTARTYEAR=-125000

export PISM_EXPERIMENT
export PISM_TITLE
export PISM_INITIALIZATION_COUPLER
export PISM_FTT_COUPLER
export PISM_FORCING_COUPLER
export PISM_REFERENCE_COUPLER
export PISM_PALEOSTARTYEAR
export PISM_FORCING_BCFILE
sh run.sh $1 $2