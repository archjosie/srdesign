#ifndef GAUSSIANBEAM_H
#define GAUSSIANBEAM_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <complex>

#include <TCanvas.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TMath.h>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TPad.h>
#include <TAxis.h>
#include <TF2.h>
#include <TH2.h>
#include <TCutG.h>
#include <TRandom.h>
#include <TColor.h>
#include <TExec.h>

#include<fftw3.h>

using namespace std;

class GaussianBeam {

private:
	static const double PI;
    
	double w0, lambda, k, zR, xMin, xMax, xInt, yMin, yMax, yInt, zMin, zMax, zInt;
	unsigned int m, n;
    vector<vector<vector<double> > > ReEField;
    vector<vector<vector<double> > > ImEField;
    vector<double> xVals, yVals, zVals;
    
	double calculateWaist(double z);
	double calculateRadCurv(double z);
	double calculateRayleigh();
	double calculateGouy(double z);
	double calculateHermite(double x, unsigned int m);
	
public:
	GaussianBeam();
	GaussianBeam(double w0, double lambda, unsigned int m, unsigned int n);
	void calculateGaussData();
	void fourierTran(vector<vector<vector<double> > > realPart, vector<vector<vector<double> > > imagPart, vector<double> xVals, vector<double> yVals, ofstream &fout);
	static double distance(double x, double y);
	void setRealE(vector<vector<vector<double> > > realE);
	void setImE(vector<vector<vector<double> > > imE);
	void setRealFour(vector<vector<vector<double> > > realFour);
	void setImFour(vector<vector<vector<double> > > imFour);
    
    void rootGraph(int argc, char** argv, unsigned int zsteps);
};

#endif
