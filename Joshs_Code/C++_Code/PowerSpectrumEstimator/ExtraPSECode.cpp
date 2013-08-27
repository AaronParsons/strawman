
void QMatrixReconstruction(NoiseCovariance& N, string outfilename){
	cout << endl << "Now reconstructing Q from basis vectors..." << endl;
	//vector< vector< vector<double> > > reconstructed(kParaBins*kPerpBins, vector< vector<double> >(xBins*yBins*fBins, vector<double>(xBins*yBins*fBins,0)));
	vector< vector<double> > reconstructed(kParaBins*kPerpBins, vector<double>(fBins, 0));
	for (int k1 = 0; k1 < fBins; k1++){
		cout << " " << floor(100.0 * k1/fBins) << "% done.\r" << std::flush;
		for (int i1 = 0; i1 < 1; i1++){
			for (int j1 = 0; j1 < 1; j1++){
/* 				for (int i2 = 0; i2 < xBins; i2++){
					for (int j2 = 0; j2 < yBins; j2++){
						for (int k2 = 0; k2 < fBins; k2++){ */
				for (int i2 = 0; i2 < 1; i2++){
					for (int j2 = 0; j2 < 1; j2++){
						for (int k2 = 0; k2 < 1; k2++){
							CVector unit1(s);
							CVector unit2(s);
							int here1 = i1*yBins*fBins + j1*fBins + k1;
							int here2 = i2*yBins*fBins + j2*fBins + k2;
							unit1.real[here1] = 1.0;	
							unit2.real[here2] = 1.0;
							CVector P = bandPowerSpectrum(unit1, unit2, N, true);
							for (int n = 0; n < kParaBins*kPerpBins; n++) reconstructed[n][here1]/*[here2]*/ = 2*P.real[n];
						}
					}
				}
			}
		}
	}
	cout << "Done.                  " << endl;	
	
	ofstream outfile;
	//outfile.precision(30);
	outfile.open(outfilename.c_str(), ios::trunc);	
	for (int n = 0; n < kParaBins*kPerpBins; n++){
		for (int i = 0; i < fBins; i++){
			//for (int j = 0; j < xBins*yBins*fBins; j++){
				outfile << reconstructed[n][i]/*[j]*/ << " ";
			//}
			//outfile << endl;
		}
		outfile << endl << endl;
	}
	outfile.close();
}
