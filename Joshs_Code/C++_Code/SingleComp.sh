#!/bin/sh

#kill all running powerSpectrumEstimators
pkill -f "powerSpectrumEstimator"

#change over all specifications files
cd Specifications
rm powerSpectrumSpecs.txt
cp powerSpectrumSpecs_$1.txt powerSpectrumSpecs.txt
rm NSpecs.txt
cp NSpecs_$1.txt NSpecs.txt
cd ..

#Calculate the noise covariance matrix
cd ./CovarianceMatrices/N/
make
./N
cd ../..

#Prepare for estimating power spectra
rm -R -f QuadraticEstimators_$1
mkdir QuadraticEstimators_$1
for ((  i = 0 ;  i <= ($2-1);  i++  ))
do
	mkdir QuadraticEstimators_$1/$i
done
cp -R ./Specifications ./QuadraticEstimators_$1/Specifications
cp ./cubeParameters.txt ./QuadraticEstimators_$1/cubeParameters.txt
cd QuadraticEstimators_$1
mkdir Spherical
mkdir Logs
for ((  i = 0 ;  i <= ($2-1);  i++  ))
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
for ((  i = 0 ;  i <= ($2-1);  i++  ))
do
	echo "Starting the instance $i..."
	nohup nice -n 15 ./powerSpectrumEstimator $i $2 >& ../QuadraticEstimators_$1/Logs/log_$i.txt &
	sleep .2
done
