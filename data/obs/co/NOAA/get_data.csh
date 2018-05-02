#!/bin/csh -f
# Run this script to get the NOAA CO data
rm -rf ./month ; mkdir ./month
wget -r --no-parent -A 'co_*_month.txt' ftp://aftp.cmdl.noaa.gov/data/trace_gases/co/flask/surface/
mv aftp.cmdl.noaa.gov/data/trace_gases/co/flask/surface/co_*_month.txt ./month/.
rm -rf aftp.cmdl.noaa.gov
exit(0)
