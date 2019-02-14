#!/bin/sh

# script to download all the OCDB files of a given year
# for some defined types.
# The intent being to have a local copy of all the files
# required to do a scan of the OCDB, e.g. to compute
# the integrated luminosity, or assess the HV issues
# of the year, etc...

year=${1:-2018}

#grp="GRP/CTP/Config GRP/GRP/Data GRP/CTP/Scalers GRP/GRP/LHCData"
grp=""

# muon="MUON/Calib/MappingData MUON/Calib/MappingRunData MUON/Calib/RecoParam MUON/Calib/Pedestals MUON/Calib/Config MUON/Calib/HV MUON/Calib/LV MUON/Calib/BPEVO MUON/Calib/OccupancyMap MUON/Calib/RejectList"
muon="MUON/Calib/OccupancyMap"

match="/"
repl="_"

# generate the lists first
for t in $(echo $grp $muon)
do
    bpath=/alice/data/$year/OCDB/$t
    list=${bpath//$match/$repl}.txt
    if test -e $list; then
        echo "File list $list already exist. Not replacing it. Do it by hand if need be."
    else
        echo "Making file list for $t -> $list"
        alien_find $bpath "*.root" | grep alice | sed s_/alice_alien:///alice_g >  $list
    fi
done

# now loop over the lists to download files
for t in $(ls _alice_data_${year}_*.txt)
do
    echo "Downloading files for $t"
    root -b <<EOF
gSystem->Load("$HOME/github.com/aphecetche/aafu/libmyaf.so");
VAF::CopyFromRemote("$t");
.q
EOF
done
