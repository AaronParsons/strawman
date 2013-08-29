#!/bin/sh



#cp $dataDIR
#set currentDIR = `pwd | awk '{print $1}'`
#set dataDIR = {$currentDIR}'/rrkltnlkertner'

echo "How many instances should I run?"
read instances


#kill all running powerSpectrumEstimators
pkill -f "powerSpectrumEstimator"

#Calculate the resolved point source covariance matrix
cd ./CovarianceMatrices/R/
#echo "Re-roll resolved point sources? (y/n)"
#read reroll
#if [ "$reroll" == "y" ]
#then
#    make
#    ./R
#fi

#Calculate the noise covariance matrix
cd ../N
make
./N

#Calculate the unresolved point sources covariance matrix
#cd ../U/
#make
#./U
#
##Calculate the galactic synchrotron covariance matrix
#cd ../G/
#make
#./G

cd ../..

#Prepare for estimating power spectra
rm -R -f QuadraticEstimators
mkdir QuadraticEstimators
for ((  i = 0 ;  i <= ($instances-1);  i++  ))
do
	mkdir QuadraticEstimators/$i
done
cp -R ./Specifications ./QuadraticEstimators/Specifications
cp ./cubeParameters.txt ./QuadraticEstimators/cubeParameters.txt
cd QuadraticEstimators
mkdir Spherical
mkdir Logs
for ((  i = 0 ;  i <= ($instances-1);  i++  ))
do
	mkdir Spherical/$i
done
rm -f iterationCounter.txt
cp ../PowerSpectrumEstimator/copiedCounter.txt ./iterationCounter.txt
cd ..



#Estimate power spectrum and Fisher matrix
cd ./PowerSpectrumEstimator
make
echo "Now Estimating the Power Spectrum and Fisher Matrix..."
for ((  i = 0 ;  i <= ($instances-1);  i++  ))
do
	echo "Starting the instance $i..."
	nohup nice -n 5 ./powerSpectrumEstimator $i >& ../QuadraticEstimators/Logs/log_$i.txt &
done

cd ..
while true
do
	sleep 1
	more QuadraticEstimators/iterationCounter.txt
done
