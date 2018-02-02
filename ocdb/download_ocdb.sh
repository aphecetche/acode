#!/bin/sh

# script to download all the OCDB files of a given year
# for some defined types.
# The intent being to have a local copy of all the files
# required to do a scan of the OCDB, e.g. to compute 
# the integrated luminosity, or assess the HV issues
# of the year, etc...

year=2017

#for type in GRP/CTP/Config GRP/GRP/Data GRP/CTP/Scalers GRP/GRP/LHCData 
for type in MUON/Calib/Pedestals MUON/Calib/Config MUON/Calib/HV MUON/Calib/LV MUON/Calib/BPEVO MUON/Calib/OccupancyMap
do
    path=/alice/data/$year/OCDB/$type
    list=${path////_}.txt
    echo "Making file list for $type -> $list"
    alien_find $path "*.root" | grep alice | sed s_/alice_alien:///alice_g >  $list
    echo "Downloading files for $type"
    root -b <<EOF
gSystem->Load("$HOME/github.com/aphecetche/aafu/libmyaf.so");
VAF::CopyFromRemote("$list");
.q
EOF
done
