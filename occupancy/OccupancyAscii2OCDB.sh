
for f in $(ls occupancy.*) 
do
  run=$(echo $f | cut -c 11-20)
  echo $run
  aliroot -q -b "OccupancyAscii2OCDB.C($run)"
done

