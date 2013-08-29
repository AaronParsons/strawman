#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include "../../CommonClasses/Specs.h"
#include <time.h>
#include "../../CommonClasses/CVector.h"

using namespace std;

//Constants and global variables
const double pi = 3.1415926535897932384626433832795;
const double c = 299792000; //m/s
const double H0 = 70900; //m/s/Mpc
const double OmegaM = .27;
const double OmegaL = .73;
const double f21cm = 1420.41; //MHz
double deltaRedshift = .00001;
int comovingDistIntSteps = 10000;
int rotationSynthesisTimeSteps = 101;
Specs *s;
double antennaDiameter, systemTemperatureConstant, systemTemperatureCoeff, systemTemperautreExp, observationDays, observationHoursPerDay, centralHourAngle, latitude, FWHMinDegrees, FWHMreferenceFrequency;
double effectiveArea, totalObsTimeInSec, fullFieldOfView, fLength, xyLength, fStart;
bool useEmpiricalSigma, useArrayFromFile;
string maskDataCubeFilename, arrayFilename;

int xBins, yBins, fBins, nElements, nAntennas;

void loadSpecs(string cubeParametersFilename, string NSpecsFilename){
	
	//Load in data cube parameters
	fstream infile(cubeParametersFilename.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> xBins;
	infile >> dummy >> yBins;
	infile >> dummy >> fBins;
	infile >> dummy >> xyLength;
	infile >> dummy >> fLength;
	infile >> dummy >> fStart; 
	infile.close();
	
	//Load relevant specifications into Specs object
	s->xBins = xBins;
	s->yBins = yBins;
	s->fBins = fBins;
	s->fLength = fLength;
	s->fStart = fStart;
	nElements = xBins*yBins*fBins;
	
	//Load in noise specs
	infile.open(NSpecsFilename.c_str(),fstream::in);
	infile >> dummy >> antennaDiameter;
	infile >> dummy >> systemTemperatureConstant;
	infile >> dummy >> systemTemperatureCoeff;
	infile >> dummy >> systemTemperautreExp;
	infile >> dummy >> observationDays;
	infile >> dummy >> observationHoursPerDay;
	infile >> dummy >> arrayFilename;
	cout << arrayFilename << endl;
	infile >> dummy >> centralHourAngle;
	infile >> dummy >> latitude;
	infile >> dummy >> FWHMinDegrees;
	infile >> dummy >> FWHMreferenceFrequency;
	infile.close();	
}


void printObsTimesToFile(vector< vector<double> >& observationTimes){
	ofstream outfile;
	string ObservationTimesFilename = "obsTimes.dat";
	outfile.open(ObservationTimesFilename.c_str(), ios::trunc);	
	for (int j = 0; j < yBins; j++){
		for (int i = 0; i < xBins; i++) outfile << observationTimes[i][j] << " ";
		outfile << endl;
	}
	outfile.close();		
}

struct coord {
	double x,y;
	bool operator<(const coord& a) const { 
		if (x*x+y*y < a.x*a.x+a.y*a.y) return true;
		if (x*x+y*y == a.x*a.x+a.y*a.y && atan2(x,y) < atan2(a.x,a.y)) return true;
		else return false;
	}
};

vector<coord> layoutFromFile(){
	fstream infile(arrayFilename.c_str(),fstream::in);
	infile >> nAntennas;
	vector<coord> antennae;
	for (int n = 0; n < nAntennas; n++){
		coord current;
		infile >> current.x;
		infile >> current.y;
		antennae.push_back(current);
	}
	return antennae;
}

map<coord,int> calculateBaselines(vector<coord>& antennae){
	cout << "Now calculating baselines." << endl;
	int percent = 0;
	map<coord,int> baselines;
	for (int i = 0; i < nAntennas; i++){
		for (int j = 0; j < nAntennas; j++){
			int n = i*nAntennas + j;
			if (100.0*n/nAntennas/nAntennas > percent){
				percent = int(ceil(100.0*n/nAntennas/nAntennas));
				cout << " " << percent << "% complete. \r" << std::flush;
			}
			coord difference;
			difference.x = antennae[j].x - antennae[i].x;
			difference.y = antennae[j].y - antennae[i].y;
			//cout << difference.x << " " << difference.y << endl;
			double baselineLength = sqrt(pow(difference.x,2) + pow(difference.y,2));
			if (!(difference.x == 0 && difference.y == 0)) baselines[difference]++;
		}
	}
	return baselines;
}

double uv2kConversion(double f){
	double comovingDist = 0;
	double z = f21cm/f - 1;
	for (int i = 0; i < comovingDistIntSteps; i++){
		double zLeft = z*i/comovingDistIntSteps;
		double zRight = z*(i+1)/comovingDistIntSteps;
		comovingDist += (1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*z/comovingDistIntSteps/6;
	}
	comovingDist *= c/H0;
	return 2*pi/comovingDist; //TODO: decide if there needs to be an extra 2pi factor here
}

vector< vector<double> > rotationSynthesis(map<coord,int>& baselines, double f, double uv2k, int slice){
	
	//FOR DRIFT SCAN INSTRUMENT

	double startHA = centralHourAngle*2*pi/24 - fullFieldOfView*2*pi/360/2;
	double endHA = centralHourAngle*2*pi/24 + fullFieldOfView*2*pi/360/2;
	double deltaHA = (endHA - startHA)/rotationSynthesisTimeSteps;

	vector< vector<double> > observationTimes(xBins, vector<double>(yBins,0));
	map<coord,int>::iterator it;
	int percent = 0;
	int mapSize = baselines.size();
	int n = 0;
	cout << endl << "Now performing rotation synthesis for frequency slice " << slice << "." << endl;
	for (it=baselines.begin() ; it != baselines.end(); it++ ){
		if (100.0*n/mapSize > percent) {
			percent = int(ceil(100.0*n/mapSize));
			cout << " " << percent << "% complete. \r" << std::flush;
		}
		n = n+1;
		double Xlambda = (*it).first.y * -sin(latitude*2*pi/360) / (c/f/1000000);
		double Ylambda = (*it).first.x / (c/f/1000000);
		double Zlambda = (*it).first.y * cos(latitude*2*pi/360) / (c/f/1000000);
		for (double ha = startHA + deltaHA/2; ha <= (endHA - deltaHA/2); ha += deltaHA){
			double kx = (Xlambda*sin(ha) + Ylambda*cos(ha))*uv2k;
			double ky = ((-Xlambda*cos(ha) + Ylambda*sin(ha)) * sin(latitude*2*pi/360) + Zlambda*cos(latitude*2*pi/360))*uv2k;
			double deltaK = 2*pi/xyLength;

			int uBinAbove = int(ceil(kx/deltaK) + yBins/2);
            int uBinBelow = int(floor(kx/deltaK) + yBins/2);
            int vBinAbove = int(ceil(ky/deltaK) + yBins/2);
            int vBinBelow = int(floor(ky/deltaK) + yBins/2);

            double uAboveRatio = 1-(uBinAbove - (kx/deltaK + xBins/2));
            double vAboveRatio = 1-(vBinAbove - (ky/deltaK + yBins/2));

			if ((uBinAbove >= 0) && (uBinAbove < xBins) && (vBinAbove >= 0) && (vBinAbove < yBins)) observationTimes[uBinAbove][vBinAbove] += uAboveRatio*vAboveRatio * (*it).second*totalObsTimeInSec/rotationSynthesisTimeSteps;
            if ((uBinAbove >= 0) && (uBinAbove < xBins) && (vBinBelow >= 0) && (vBinBelow < yBins)) observationTimes[uBinAbove][vBinBelow] += uAboveRatio*(1-vAboveRatio) * (*it).second*totalObsTimeInSec/rotationSynthesisTimeSteps;
            if ((uBinBelow >= 0) && (uBinBelow < xBins) && (vBinAbove >= 0) && (vBinAbove < yBins)) observationTimes[uBinBelow][vBinAbove] += (1-uAboveRatio)*vAboveRatio * (*it).second*totalObsTimeInSec/rotationSynthesisTimeSteps;
            if ((uBinBelow >= 0) && (uBinBelow < xBins) && (vBinBelow >= 0) && (vBinBelow < yBins)) observationTimes[uBinBelow][vBinBelow] += (1-uAboveRatio)*(1-vAboveRatio) * (*it).second*totalObsTimeInSec/rotationSynthesisTimeSteps;
		}
	}
	return observationTimes;
}

vector< vector< vector<double> > > calculateObservationTimesFromArray(vector<double>& freqs){
	vector<coord> antennae = layoutFromFile();
	map<coord,int> baselines = calculateBaselines(antennae);

	vector< vector< vector<double> > > observationTimes;
	for (int k = 0; k < fBins; k++){
		double uv2k = uv2kConversion(freqs[k]);
		vector< vector<double> > observationTimeSlice = rotationSynthesis(baselines, freqs[k], uv2k, k);
		observationTimes.push_back(observationTimeSlice);
	}
	
	double tMax = 0;
	for (int k = 0; k < fBins; k++){
		for (int j = 0; j < yBins; j++){
			for (int i = 0; i < xBins; i++){
				if (observationTimes[k][i][j] > tMax) tMax = observationTimes[k][i][j];
			}
		}
	}
/*	for (int k = 0; k < fBins; k++){
		for (int j = 0; j < yBins; j++){
			for (int i = 0; i < xBins; i++){
				if (observationTimes[k][i][j] < tMax * leastObservedCutoff){
					for (int k2 = 0; k2< fBins; k2++) observationTimes[k2][i][j] = 0;
				}
			}
		}
	}*/
	

	printObsTimesToFile(observationTimes[fBins/2]);
	return observationTimes;
}

vector<double> listFrequencies(){
	double comovingDist = 0;
	double zRight = f21cm/fStart - 1;
	double zLeft = zRight + deltaRedshift;
	while (true){
		comovingDist -= c/H0*((1.0/sqrt(OmegaM*pow(1+zRight,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zLeft,3)+OmegaL))*deltaRedshift/6);
		if (comovingDist <= -fLength) break;
		zRight = zLeft;
		zLeft = zRight - deltaRedshift;
	}
	double fL = f21cm/(zLeft + 1);
	double fEnd = fL + (fStart - fL)/fBins;
	vector<double> freqs(fBins,0);
	double deltaF = (fEnd - fStart)/(fBins - 1);
	for (int i = 0; i < fBins; i++){
		freqs[i] = fStart + (fBins - 1 - i)*deltaF;
	} 
	return freqs;
}

double calculateOmegaPix(vector<double>& freqs){	
	double f = fStart + (freqs[fBins -1] - fStart)/2;		
	double comovingDist = 0;
	double z = f21cm/f - 1;
	for (int i = 0; i < comovingDistIntSteps; i++){
		double zLeft = z*i/comovingDistIntSteps;
		double zRight = z*(i+1)/comovingDistIntSteps;
		comovingDist += (1.0/sqrt(OmegaM*pow(1+zLeft,3) + OmegaL) + 4.0/sqrt(OmegaM*pow(1+(zLeft + zRight)/2,3) + OmegaL) + 1.0/sqrt(OmegaM*pow(1+zRight,3)+OmegaL))*z/comovingDistIntSteps/6;
	}
	comovingDist *= c/H0;
	double Omega = xyLength * xyLength / comovingDist / comovingDist;
	return Omega / xBins / yBins;
}

double calculateOmegaPrime(vector<double>& freqs){	
	double sigma = fullFieldOfView * 2 * pi / 360 / 2.3548201;
	double OmegaP = 2*pi*pow(sigma,2)*pow(erf(pi/sqrt(2)/sigma),2);
	double OmegaPP = pi * pow(sigma,2)* pow(erf(pi/sigma),2);
	return pow(OmegaP,2)/OmegaPP;
}

/*double calculateTMax(vector< vector< vector<double> > >& observationTimes){
	double tMax = 0;
	for (int k = 0; k < fBins; k++){
		for (int j = 0; j < yBins; j++){
			for (int i = 0; i < xBins; i++){
				if (observationTimes[k][i][j] > tMax) tMax = observationTimes[k][i][j];
			}
		}
	}
	return tMax;
}*/


void updateN(int n, double deltaNu, double OmegaPix, double OmegaPrime, double f, vector< vector<double> >& observationTimes, vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	double Tsys = systemTemperatureConstant + systemTemperatureCoeff * pow(c/1e6/f,systemTemperautreExp);
	cout << "For slice " << n << ", Tsys = " << Tsys << endl;
	double prefactor = pow(Tsys,2) * OmegaPrime / 2 / (deltaNu*1e6) / OmegaPix;
	for (int i = 0; i < xBins; i++){
		for (int j = 0; j < yBins; j++){
			if (observationTimes[i][j] == 0 || i == 0 || j == 0) noiseCovarianceMatrixDiagonal[i][j][n] = -1;
			else noiseCovarianceMatrixDiagonal[i][j][n] = prefactor / observationTimes[i][j];
		}
	}
}

void printNtoFile(vector< vector < vector<double> > >& noiseCovarianceMatrixDiagonal){
	string NOutputFilename = "N.dat";
	ofstream outfile;
	outfile.precision(30);
	outfile.open(NOutputFilename.c_str(), ios::trunc);	
	for (int u = 0; u < xBins; u++){
		for (int v = 0; v < yBins; v++){
			for (int k = 0; k < fBins; k++){
				outfile << noiseCovarianceMatrixDiagonal[u][v][k] << endl;
			}
		}
	}
	outfile.close();
}

int main(){
	//Load in parameters and data
	cout << endl << "Now calculating the noise covariance matrix..." << endl;
	s = new Specs();
	loadSpecs("../../cubeParameters.txt","../../Specifications/NSpecs.txt");

	//Calculate wavelengths and geometric factors
	vector<double> freqs = listFrequencies();
	double deltaF = freqs[0] - freqs[1];
	double bandwidth = fBins * deltaF; 
	fullFieldOfView = FWHMinDegrees * FWHMreferenceFrequency / freqs[fBins/2];
	double omegaPix = calculateOmegaPix(freqs);
	double OmegaPrime = calculateOmegaPrime(freqs); //Note: the f depenedence of omegaPix and OmegaPrime are assumed to cancel
	
	//Calculate observation times
	totalObsTimeInSec = 3600.0 * observationDays * observationHoursPerDay * fullFieldOfView / 360;
	cout << "This " << fullFieldOfView << " degree by " << fullFieldOfView << " degree field is observed for " << totalObsTimeInSec/3600 << " hours in " << observationDays << " days of operation." << endl;
	vector< vector< vector<double> > > observationTimes = calculateObservationTimesFromArray(freqs);
	
	//Calculate N and save
	vector< vector < vector<double> > > noiseCovarianceMatrixDiagonal(xBins, vector< vector<double> >(yBins, vector<double>(fBins,0)));
	for (int n = 0; n < fBins; n++) updateN(n,deltaF,omegaPix,OmegaPrime,freqs[n],observationTimes[n],noiseCovarianceMatrixDiagonal);
	printNtoFile(noiseCovarianceMatrixDiagonal);
	cout << "Noise covariance calculation complete." << endl << endl;
	return 0;
}
