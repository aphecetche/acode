#!/bin/sh
#
#
# stress test to execute from within the root/src/test directory
#
#

rm -rf dir.stress.*
rm -rf result.*

# do each test nruns times to get "some" statistical significance ?
nruns=1 
withio=1
ncores=$(getconf _NPROCESSORS_ONLN)

for ((run=0;run<$nruns;++run))
do

	a=$((run+1))
	printf "Run %d/%d : " $a $nruns

	for ((n=1;n<=$ncores;++n))
	do

		for ((i=0;i<$n;++i)) do
			mkdir -p dir.stress.$i
			cd dir.stress.$i
			../rootstress $n $withio ../files > stress.$n.log &
			cd ..
		done

		FAIL=0

		for job in `jobs -p`
		do
    		wait $job || let "FAIL+=1"
		done

		if [ "$FAIL" == "0" ];
		then
			printf "ALL OK for %d processes - " $n
		else
			echo "FAILURE ! ($FAIL)"
		fi

	done

	printf "\n"

	cat dir.stress.0/stress.*.log | cut -f 1 -d '*' > result.$run

done

root -b -q plot.C++\($nruns\)