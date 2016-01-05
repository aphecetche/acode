#!/bin/sh

#usage:
#$1 = base directory where to look for the logs (e.g. /alice/data/2010/LHC10h/000137161/ESDs/pass1)
#  or *.txt file containing the list of chunks to check (e.g. /alice/data/2011/LHC11h/000168207/ESDs/pass2_muon/11000168207078.29)
#$2 = output directory
#$3 = if present, assume the input is a text file with xrootd accessible log archives (e.g. root://nansafmaster2.in2p3.fr//..../log_archive.zip)

if [ $# = 0 ]; then
  echo "give the base directory where to look for the logs in arguments"
  exit 3
fi

nonalien=0

if [ $# = 3 ]; then
  nonalien=1
fi

inDir=""

if [ $nonalien -eq 0 ]; then
  inDir="alien:/"
  source /tmp/gclient_env_$UID
fi

outdir="CheckLogs"
if [ $# = 2 ]; then
  outdir=$2
fi
if [[ -d $outdir ]]; then
  echo "a directory containing checks ($outdir/) already exist. Please move or remove it"
  exit 3
fi
mkdir $outdir

inFile=`pwd`/$1
if [[ $1 != *.txt ]]; then
  inDir="alien://$1"
  inFile="__chunk.txt__"
  alien_ls -F $1/ | grep "/$" | sed 's:/$::' > $outdir/$inFile
fi

cd $outdir

num=0
nGoodEv=0
nBadEvAccess=0
nBadEvBranch=0
nBadEvRawData=0
nBadEvReadout=0
nBadEvCTP=0
nBadEv=0
nChunksRecover=0

for dir in `cat $inFile`; do

  let num+=1

# if [ $num -lt 2990 ]
# then
#   continue
# fi

  echo "checking file number $num (in $inDir/$dir)"

  if [ $nonalien -eq 0 ]; then
    alien_cp $inDir/${dir}/rec.log tmp.txt
  else
    xrdcp $dir .
    unzip -p log_archive.zip > tmp.txt
    rm log_archive.zip
  fi

  if [[ ! -e tmp.txt ]]; then
    continue
  fi

  cat tmp.txt | grep "E-AliMUON" >> error.txt
  cat tmp.txt | grep "W-AliMUON" >> warning.txt
  cat tmp.txt | grep "According to mask 400be9b (human readable form below)" >> badPads.txt

  eval $(awk -v nGoodEvInFile=0 -v badEvAccess=0 -v badEvAccessNum="" -v badEvBranch=0 -v badEvBranchNum="" -v badEvRawData=0 -v badEvRawDataNum="" -v badEvReadoutNum="" -v badEvCTP=0 -v badEvCTPNum="" -v badEvNum="" '{
    if (index($0,"E-TAlienFile::ReadBuffer: The remote file is not open") > 0) badEvAccess=1;
    else if (index($0,"E-TBranchElement::GetBasket") > 0 || index($0,"E-TBranchRef::GetBasket") > 0) badEvBranch=1;
    else if (index($0,"of type LowMultiplicity ===") > 0 || index($0,"of type HighMultiplicity ===") > 0) {
      nGoodEvInFile++;
      while ((getline line) > 0 && index(line,"=== End Event") == 0) {
        if (index(line,"shower event") > 0) print line >> "triggerPb.txt";
        else if (index(line,"Too many track candidates") > 0) print line >> "trackerPb.txt";
        else if (index(line,"Too many local maxima") > 0) print line >> "clusterPb.txt";
        else if (index(line,"E-TAlienFile::ReadBuffer: The remote file is not open") > 0) badEvAccess=1;
        else if (index(line,"E-TBranchElement::GetBasket") > 0 || index(line,"E-TBranchRef::GetBasket") > 0) badEvBranch=1;
        else if (index(line,"raw data size found in the header is wrong") > 0) badEvRawData=1;
        else if (index(line,"No valid CTP (trigger) DDL raw data is found") > 0) badEvCTP=1;
      }
      if (badEvAccess == 1) badEvAccessNum=badEvAccessNum","$5;
      if (badEvBranch == 1) badEvBranchNum=badEvBranchNum","$5;
      if (badEvRawData == 1) badEvRawDataNum=badEvRawDataNum","$5;
      if (badEvAccess == 1 || badEvBranch == 1 || badEvRawData == 1) badEvReadoutNum=badEvReadoutNum","$5;
      if (badEvCTP == 1) badEvCTPNum=badEvCTPNum","$5;
      if (badEvAccess == 1 || badEvBranch == 1 || badEvRawData == 1 || badEvCTP == 1) badEvNum=badEvNum","$5;
      badEvAccess=0;
      badEvBranch=0;
      badEvRawData=0;
      badEvCTP=0;
    } else if (index($0,"=== End Event") > 0) {
      badEvAccess=0;
      badEvBranch=0;
    }
  } END {
    print "nGoodEvCurrent="nGoodEvInFile";";
    print "badEvAccessNumCurrent=\""badEvAccessNum"\";";
    print "badEvBranchNumCurrent=\""badEvBranchNum"\";";
    print "badEvRawDataNumCurrent=\""badEvRawDataNum"\";";
    print "badEvReadoutNumCurrent=\""badEvReadoutNum"\";";
    print "badEvCTPNumCurrent=\""badEvCTPNum"\";";
    print "badEvNumCurrent=\""badEvNum"\";";
  }' tmp.txt)

  let "nGoodEv += nGoodEvCurrent"
  echo $nGoodEv > nEvents.txt

  if [[ $badEvAccessNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvAccessNumCurrent/,/}" >> accessPb.txt
    let "nBadEvAccess += `echo ${badEvAccessNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEvAccess > nBadEventsAccess.txt
  fi

  if [[ $badEvBranchNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvBranchNumCurrent/,/}" >> branchPb.txt
    let "nBadEvBranch += `echo ${badEvBranchNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEvBranch > nBadEventsBranch.txt
  fi

  if [[ $badEvRawDataNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvRawDataNumCurrent/,/}" >> rawDataPb.txt
    let "nBadEvRawData += `echo ${badEvRawDataNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEvRawData > nBadEventsRawData.txt
  fi

  if [[ $badEvReadoutNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvReadoutNumCurrent/,/}" >> readoutPb.txt
    let "nBadEvReadout += `echo ${badEvReadoutNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEvReadout > nBadEventsReadout.txt
  fi

  if [[ $badEvCTPNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvCTPNumCurrent/,/}" >> CTPPb.txt
    let "nBadEvCTP += `echo ${badEvCTPNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEvCTP > nBadEventsCTP.txt
  fi

  if [[ $badEvNumCurrent ]]; then
    echo "$inDir/$dir/AliESDs.root@${badEvNumCurrent/,/}" >> eventPb.txt
    let "nBadEv += `echo ${badEvNumCurrent//[^,]/} | grep -o ',' | wc -l`"
    echo $nBadEv > nBadEvents.txt
  fi

  nChunksRecoverCurrent=`grep -c 'recovered key TTree:RAW' tmp.txt`
  if [[ $nChunksRecoverCurrent != 0 ]]; then
    echo "$inDir/$dir/AliESDs.root" >> recoverPb.txt
    let "nChunksRecover += nChunksRecoverCurrent"
    echo $nChunksRecover > nChunksRecover.txt
  fi

  rm -f tmp.txt

done

cd ..

if [[ $1 != *.txt ]]; then
  rm -f $outdir/$inFile
fi
