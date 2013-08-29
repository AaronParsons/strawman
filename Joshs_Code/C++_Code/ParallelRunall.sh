#!/bin/sh

echo "How many instances should I per computer?"
read instances

for ((  comp = 4 ;  comp <= 7;  comp++  ))
do
	ssh eor-0$comp "cd ~/HERA; ./SingleComp.sh $comp $instances"
	sleep 10
done

while true
do
	more QuadraticEstimators_4/iterationCounter.txt
	sleep 10
done