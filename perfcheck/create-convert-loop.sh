#!/bin/sh

echo "root -b <<EOF"
echo ".L /source/QuickAccEff.C+"

for file in $(find /data/simulations/IdealSimulations4/ -name AliESDs.root)
do
    dest=${file/AliESDs/compact}
    list="$list $dest"
    echo "ConvertESD(\"$file\",\"$dest\");"
done
echo "EOF"
echo hadd -f compact.root $list

