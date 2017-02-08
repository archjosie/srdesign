#include "GaussianBeam.h"


//const static double PI = 3.14159265;
using namespace std;

GaussianBeam::GaussianBeam() {
	lambda = 1;
	w0 = 20;
	p = 0;
	l = 0;

	vector<vector<vector<double> > > local(1, vector<vector<double> >(1, vector<double>(1, 0)));
	this->ReEField = local;
	this->ImEField = local;

	this->k = 2 * 3.14159265 / lambda;
	this->zR = calculateRayleigh();
}

GaussianBeam::GaussianBeam(double w0, double lambda, unsigned int p, unsigned int l) {
    this->w0 = w0;
	this->lambda = lambda;
	this->p = p;
	this->l = l;

	vector<vector<vector<double> > > local(1, vector<vector<double> >(1, vector<double>(1, 0)));
	this->ReEField = local;
	this->ImEField = local;

	this->k = 2 * 3.14159265 / lambda;
	this->zR = calculateRayleigh();
}

double GaussianBeam::distance(double x, double y) {
	return sqrt(pow(x, 2) + pow(y, 2));
}

double GaussianBeam::calculateRadCurv(double z) {
	return (z*(1 + pow(zR / z, 2)));
}

double GaussianBeam::calculateGouy(double z) {
	return atan(z / zR)*(abs(l) + 2*p);
}

double GaussianBeam::calculateWaist(double z) {
	return w0*sqrt(1 + pow(z / zR, 2));
}

double GaussianBeam::calculateRayleigh() {
	return 3.14159265 * pow(w0, 2) / lambda;
}

double GaussianBeam::calculateHermite(double x, unsigned int m) { //Calculates the nth Hermite polynomial evaluated at x recursively
	if (m == 0) return 1;
	if (m == 1) return 2 * x; //Base cases

	return (2 * x*calculateHermite(x, m - 1) - 2 * (m - 1)*calculateHermite(x, m - 2));
}

void GaussianBeam::calculateGaussData() {
    rMax=3;
    rMin=0;
	rInt=.1;
	tMax = 2 * 3.14159265;
	tMin = 0;
	tInt = 3.14159265 / 100;

	int rRange = (rMax - rMin) / rInt + 1;
	int tRange = (tMax - tMin) / tInt + 1

	for (int i = 0; i < rRange; i++) { //Populate rVals
		double rCurr = rMin + i*rInt;
		rVals.push_back(rCurr);
	}

	for (int i = 0; i < tRange; i++) { //Populate rVals
		double tCurr = tMin + i*tInt;
		tVals.push_back(tCurr);
	}

	double theZ; //Stores lower limit on distance from focus
	theZ=.1;
	zVals.push_back(theZ);

	vector<vector<vector<double> > > ReLocal(xVals.size(), vector<vector<double> >(tVals.size(), vector<double>(zVals.size(), 0))); //Separate 3D vectors to hold real and imaginary parts of E-field

    ReEField = ReLocal;
    ImEField = ReLocal;

	for (int m = 0; m < zVals.size(); ++m) {  //Weird choice, but this streamlines the program
										  //These only vary with z, so we only need to calculate them once per loop
		double gouy = calculateGouy(zVals.at(m));
		double radCurv = calculateRadCurv(zVals.at(m));
		double spotSize = calculateWaist(zVals.at(m));

		for (int i = 0; i < rVals.size(); ++i) {
			for (int j = 0; j < tVals.size(); ++j) {

				double imArg = -(k* zVals.at(m) + k*(pow(rVals.at(i), 2) / (radCurv * 2)) + l*tVals.at(j) - gouy);
				complex<double> expArg(0.0, imArg);

				complex<double> phasorOut; //Uncomment if we only want to consider real component of field			
				phasorOut = laguerre(2*pow(rVals.at(i),2) / pow(spotSize,2), abs(l), p)*pow(sqrt(2)*rVals.at(i) / spotSize, abs(l))*(1 / spotSize)*exp(-pow(rVals.at(i) / spotSize, 2))*exp(expArg);
				double realField = real(phasorOut);
				double imagField = imag(phasorOut);

				ReEField.at(i).at(j).at(m) = realField;
				ImEField.at(i).at(j).at(m) = imagField;
			}
		}
	}
}

void GaussianBeam::fourierTran(vector<vector<vector<double> > > realPart, vector<vector<vector<double> > > imagPart, vector<double> xVals, vector<double> yVals, ofstream &fout) {
	fftw_complex *in, *out;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * realPart.size() * realPart.at(0).size());
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * realPart.size() * realPart.at(0).size());

	int k = 0;
	for (int i = 0; i < realPart.size(); ++i) {
		for (int j = 0; j < realPart.at(0).size(); ++j) {
			in[k][0] = realPart.at(i).at(j).at(0);
			in[k][1] = imagPart.at(i).at(j).at(0);
			k++;
		}
	}

	fftw_plan g = fftw_plan_dft_2d(realPart.size(), realPart.at(0).size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(g);

	vector<vector<double> > realFour(realPart.size(), vector<double>(realPart.at(0).size()));
	vector<vector<double> > imagFour(imagPart.size(), vector<double>(imagPart.at(0).size()));

	k = 0;
	for (int i = 0; i < realPart.size(); ++i) {
		for (int j = 0; j < realPart.at(0).size(); ++j) {
			realFour.at(i).at(j) = out[k][0];
			imagFour.at(i).at(j) = out[k][1];
			k++;
		}
	}

	fout << "# x " << '\t' << "y " << '\t' << "Re(E(x,y)) " << '\t' << "Im(E(x,y)) " << '\t' << "Abs(E(x,y))" << endl; //Pound sign forces gnuplot to ignore line

																																			 //Write our 3D complex vector to file
	for (int i = 0; i < xVals.size(); ++i) {
		for (int j = 0; j < yVals.size(); ++j) {
				fout << xVals.at(i) << '\t' << yVals.at(j) << '\t' << realFour.at(i).at(j) << '\t' << imagFour.at(i).at(j)<< '\t' << distance(realFour.at(i).at(j), imagFour.at(i).at(j)) << endl;
			}
		}

	cout << "Data output complete! Thank you for your business!" << endl;
}

vector<vector<vector<double> > > GaussianBeam::getRealE()
{
	return ReEField;
}

vector<vector<vector<double> > > GaussianBeam::getImE()
{
	return ImEField;
}

vector<double> GaussianBeam::getXVals()
{
	return xVals;
}

vector<double> GaussianBeam::getYVals()
{
	return yVals;
}

double GaussianBeam::getK()
{
	return k;
}

double GaussianBeam::laguerre(double x, double alpha, double k) {
	if (k == 0) return 1;
	if (k == 1) return 1 + alpha - x;

	return ((2 * k + 1 + alpha - x)*laguerre(x, alpha, k) - (k + alpha)*laguerre(x, alpha, k - 1)) / (k + 1)
}

void GaussianBeam::rootGraph(int argc, char** argv, vector<vector<vector<double> > > Field){
    //Open root graphics
    TApplication theApp("App", &argc, argv);
    gStyle->SetOptStat(0);
//    gStyle->SetPalette(82);
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    TH2F *hcontz = new TH2F("hcontz","Gaussian Beam Cross Section",40,xMin,xMax,40,yMin,yMax);
    Float_t px, py;
	for (int i = 0; i < xVals.size(); ++i) {
		for (int j = 0; j < yVals.size(); ++j) {
				hcontz->Fill(xVals.at(i),yVals.at(j),Field.at(i).at(j).at(0));
		}
	}
//    hcontz->SetTitle("Cross Section of Gaussian Beam");
    hcontz->GetXaxis()->SetTitle("x"); 
    hcontz->GetYaxis()->SetTitle("y"); 
    hcontz->SetMarkerStyle(1);
    hcontz->GetXaxis()->CenterTitle(); 
    hcontz->GetYaxis()->CenterTitle(); 
    hcontz->SetMaximum(4);
    hcontz->SetMinimum(0);    
    hcontz->Draw("CONT1z");
    // Output PDF
    c1->Print("GBplots.pdf","pdf");
    theApp.Run();
}
