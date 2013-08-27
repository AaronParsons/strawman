#!/bin/sh

echo -e "\n--------------------------"

for ((  comp = 4 ;  comp <= 7;  comp++  ))
do
	echo -n "On eor-0$comp: "
	find QuadraticEstimators_$comp/[0-9]*/pse* | wc -l
	echo -n "Processes running: "
	ssh eor-0$comp "ps -C powerSpectrumEstimator | grep -c power"
done

echo -e "--------------------------\n"