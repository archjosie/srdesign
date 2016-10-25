#ifndef GAUSSIANBEAM_H
#define GAUSSIANBEAM_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <complex>

#include<fftw3.h>

using namespace std;

class GaussianBeam {

private:
	static const double PI;

	double w0, lambda, k, zR;
	unsigned int m, n;

	double calculateWaist(double z);
	double calculateRadCurv(double z);
	double calculateRayleigh();
	double calculateGouy(double z);
	double calculateHermite(double x, unsigned int m);
	
public:
	GaussianBeam();
	GaussianBeam(double w0, double lambda, unsigned int m, unsigned int n);
	void calculateGaussData();
	void fourierTran(vector<vector<vector<double> > > realPart, vector<vector<vector<double> > > imagPart, vector<double> xVals, vector<double> yVals, ofstream & fout);
	double distance(double x, double y);
};

#endif
