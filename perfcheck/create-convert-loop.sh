#!/bin/sh

echo "root -b <<EOF"
echo ".L QuickAccEff.C+"

for file in $(find /alice/cern.ch/user/l/laphecet/simulations/idealpp13/ -name AliESDs.root)
do
    dest=${file/AliESDs/compact}
    list="$list $dest"
    echo "ConvertESD(\"$file\",\"$dest\",\"local:///alice/data/2015/OCDB\");"
done
echo "EOF"
echo hadd -f compact.root $list

