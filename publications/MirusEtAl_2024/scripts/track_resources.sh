#!/bin/bash

# Tim Mirus
# track the resources used by a given process

if [ ! $# -eq 3 ]; then
	echo "Correct Usage: ./track_resources.sh <ProcessName> <Interval> <OutputFile>"
fi

Interval=$2
File=$3
pName=$1

if [ -f $File ]; then
	echo "File $File exists; Abort"
	exit
fi
touch $File

echo "Time	RSS	CPU" > $File
counter=0
SECONDS=0
while true
do
	resources=`top -b -n 8 -d.1 | grep $pName | grep -v .sh | grep -v + | tail -n 1 | awk '{print $6 "\t" $9}'`
	if [[ $resources != "" ]]; then
		echo "$SECONDS	$resources" >> $File
	fi
	
	let counter=$counter+1
	sleep $Interval
done
