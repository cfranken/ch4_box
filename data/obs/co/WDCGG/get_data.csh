#!/bin/csh -f
# This script should download the CO data from the WDCGG page.
# However, getting the CO data requires submitting a query to the WDCGG page.
# The script *should* submit the query but sometimes it's a bit slow to archive the data.
# You can manually download the tarball from the WDCGG page here: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/download.cgi?para=CO"
# Choose the "event" data here: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/archiver.cgi?tool=gzip&archive=CO_EVENT"
# This should give you a CO_EVENT.tgz tarball.
# Untar the tarball and manually run the commands from this script.

rm -rf ./event CO_EVENT.tgz ; mkdir ./event
open "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/archiver.cgi?tool=gzip&archive=CO_EVENT"
wget "http://ds.data.jma.go.jp/gmd/wdcgg/tmp/ftptmp/CO_EVENT.tgz"
tar -xvf CO_EVENT.tgz
mv co/event/* ./event/.
rm -rf co CO_EVENT.tgz
exit(0)
