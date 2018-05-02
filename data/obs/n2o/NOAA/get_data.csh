#!/bin/csh -f
# Run this script to get the NOAA N2O data
rm -rf ./combined ; mkdir ./combined
wget "ftp://ftp.cmdl.noaa.gov/hats/n2o/combined/GMD_global_N2O.txt"
mv GMD_global_N2O.txt ./combined/.
exit(0)
