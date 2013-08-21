#!/bin/sh


#kill all running powerSpectrumEstimators
pkill -f "powerSpectrumEstimator"

#Parse the data and mask into my format and compute data cube specifications
cd ../DataParser
make
./dataParser

#Calculate the noise covariance matrix
cd ../CovarianceMatrices/N
make
./N

#Calculate the unresolved point sources covariance matrix
cd ../U/
make
./U

#Calculate the galactic synchrotron covariance matrix
cd ../G/
make
./G
cd ../..

#Prepare for estimating power spectra
rm -R -f SimulatedQhats
mkdir SimulatedQhats
for ((  i = 0 ;  i <= (0);  i++  ))
do
	mkdir SimulatedQhats/$i
done
cd SimulatedQhats
mkdir Spherical
for ((  i = 0 ;  i <= (0);  i++  ))
do
	mkdir Spherical/$i
done
cd ..
rm -f iterationCounter.txt
cp ./PowerSpectrumEstimator/copiedCounter.txt ./iterationCounter.txt

#Estimate power spectrum and Fisher matrix
cd ./PowerSpectrumEstimator
make
