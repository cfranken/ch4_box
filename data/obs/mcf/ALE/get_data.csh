#!/bin/csh -f
# Run this script to get the ALE MCF data
rm ./*-ale.mon
wget -r --no-parent -A '*-ale.mon' http://agage.eas.gatech.edu/data_archive/ale/monthly/
mv agage.eas.gatech.edu/data_archive/ale/monthly/*-ale.mon ./.
rm -rf agage.eas.gatech.edu
exit(0)
