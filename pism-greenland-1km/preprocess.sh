#!/bin/bash

# Copyright (C) 2009-2011 Andy Aschwanden

# downloads ice2sea "Present Day Greenland" master dataset
# downloads SeaRISE "Present Day Greenland" master dataset
# NetCDF file, adjusts metadata, and saves under new name, ready for PISM
# depends on wget and NCO (ncrename, ncap, ncatted, ncpdq, ncks)

set -e  # exit on error

# get 1km from PISM website
for DATANAME in  "pism_Greenland_1km_griggs_updated_smoothed.nc" "pism_Greenland_2km_griggs_updated_smoothed.nc" "pism_Greenland_2.5km_griggs_updated_smoothed.nc" "pism_Greenland_5km_griggs_updated_smoothed.nc" "pism_Greenland_10km_griggs_updated_smoothed.nc" "pism_Greenland_20km_griggs_updated_smoothed.nc"; do
    wget -nc http://pism-docs.org/download/$DATANAME
done
# generate config file
echo "  Generating config files..."
for CONFIG in "initial_config" "era_config"; do
ncgen -o ${CONFIG}.nc ${CONFIG}.cdl
done
echo "  Done generating config file."
echo


# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
SRG_VERSION=0.93
SRG_PATH=http://websrv.cs.umt.edu/isis/images/8/86/
SRG_NAME=Greenland_5km_v$SRG_VERSION.nc
wget -nc $SRG_PATH/$SRG_NAME

PISMVERSION=pism_$SRG_NAME
ncks -O $SRG_NAME $PISMVERSION  # just copies over, but preserves history and global attrs

# adjust metadata; uses NCO (http://nco.sourceforge.net/)
ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y $PISMVERSION
ncrename -O -v time,t -d time,t $PISMVERSION
ncrename -O -v usrf,usurf $PISMVERSION
ncrename -O -v presartm,artm $PISMVERSION
# convert from water equiv to ice thickness change rate; assumes ice density 910.0 kg m-3
ncap -O -s "precip=presprcp*(1000.0/910.0)" $PISMVERSION $PISMVERSION
ncatted -O -a units,precip,a,c,"m a-1" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precip,mapping,artm \
  $PISMVERSION $PISMVERSION
echo "  PISM-readable file $PISMVERSION created from $SRG_NAME"
echo "    (contains only fields used in bootstrapping ...)"

TEMPSERIES=pism_dT.nc
SLSERIES=pism_dSL.nc
echo -n "creating paleo-temperature file $TEMPSERIES from $SRG_NAME for option -atmosphere ...,delta_T ... "
ncks -O -v oisotopestimes,temp_time_series $SRG_NAME $TEMPSERIES
ncrename -O -d oisotopestimes,time -v oisotopestimes,time -v temp_time_series,delta_T $TEMPSERIES
ncpdq -O --rdr=-time $TEMPSERIES $TEMPSERIES  # reverse time dimension
ncap -O -s "time=-time" $TEMPSERIES $TEMPSERIES  # make times follow same convention as PISM
ncatted -O -a units,time,a,c,"years since 1-1-1" $TEMPSERIES
ncatted -O -a units,delta_T,m,c,"Kelvin" $TEMPSERIES
echo "done."
echo
echo -n "creating paleo-sea-level file $SLSERIES from $SRG_NAME for option -ocean ...,delta_SL ... "
ncks -O -v sealeveltimes,sealevel_time_series $SRG_NAME $SLSERIES
ncrename -O -d sealeveltimes,time -v sealeveltimes,time -v sealevel_time_series,delta_SL $SLSERIES
ncpdq -O --rdr=-time $SLSERIES $SLSERIES  # reverse time dimension
ncap -O -s "time=-time" $SLSERIES $SLSERIES  # make times follow same convention as PISM
ncatted -O -a units,time,a,c,"years since 1-1-1" $SLSERIES
echo "done."
echo
