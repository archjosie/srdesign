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
	double PI = 3.14159265;


	dimset = 301;
	xMax = 400;
	// xMax=.01;
	xMin = -xMax;
	xInt = 2 * xMax / (dimset - 1);
	yMax = xMax;
	yMin = -yMax;
	yInt = xInt;

	int xRange = (xMax - xMin) / xInt + 1;
	int yRange = (yMax - yMin) / yInt + 1;

	for (int i = 0; i < xRange; i++) { //Populate xVals
		double xCurr = xMin + i*xInt;
		xVals.push_back(xCurr);
		// cout << "New x: " << xCurr << endl;
	}

	for (int i = 0; i < yRange; i++) { //Populate rVals
		double yCurr = yMin + i*yInt;
		yVals.push_back(yCurr);
		// cout << "New y: " << yCurr << endl;
	}

	double theZ; //Stores lower limit on distance from focus
	theZ = .1;
	zVals.push_back(theZ);

	vector<vector<vector<double> > > ReLocal(xVals.size(), vector<vector<double> >(yVals.size(), vector<double>(zVals.size(), 0))); //Separate 3D vectors to hold real and imaginary parts of E-field

	ReEField = ReLocal;
	ImEField = ReLocal;

	/*for (int m = 0; m < zVals.size(); ++m) {  //Weird choice, but this streamlines the program
	//These only vary with z, so we only need to calculate them once per loop
	double gouy = calculateGouy(zVals.at(m));
	double radCurv = calculateRadCurv(zVals.at(m));
	double spotSize = calculateWaist(zVals.at(m));

	for (int i = 0; i < xVals.size(); ++i) {
	for (int j = 0; j < yVals.size(); ++j) {

	double r = distance(xVals.at(i), yVals.at(j));
	double t;
	if (abs(yVals.at(j)) < PI*1e-10) t = PI - sign(xVals.at(i))*PI / 2;
	else t = atan(xVals.at(i) / yVals.at(j));

	double imArg = -(k* zVals.at(m) + k*(pow(r, 2) / (radCurv * 2)) + l*t - gouy);
	complex<double> expArg(0.0, imArg);
	//cout << expArg << endl;

	complex<double> phasorOut; //Uncomment if we only want to consider real component of field
	phasorOut = laguerre(2*pow(r,2) / pow(spotSize,2), abs(l), p)*pow(sqrt(2)*r / spotSize, abs(l))*(1 / spotSize)*exp(-pow(r / spotSize, 2))*exp(expArg);
	//cout << phasorOut << endl;
	double realField = real(phasorOut);
	double imagField = imag(phasorOut);

	ReEField.at(i).at(j).at(m) = realField;
	ImEField.at(i).at(j).at(m) = imagField;
	}
	}
	}*/

	double omega = 20 / k;

	for (int i = 0; i < xVals.size(); ++i) {
		for (int j = 0; j < yVals.size(); ++j) {
			double r = distance(xVals.at(i), yVals.at(j));

			double reArg = -pow(k, 2)*pow(r, 2) / (pow(omega*k, 2));
			double imArg = 0;
			if (abs(yVals.at(j)) < PI*1e-10) imArg = -l*(PI - sign(xVals.at(i))*PI / 2);
			else imArg = -(l*atan(xVals.at(i) / yVals.at(j)));
			//cout << reArg << endl;
			
			complex<double> phasorOut; //Uncomment if we only want to consider real component of field			
			phasorOut = pow(2, abs(l) / 2)/ (omega * k) * exp(reArg+100) * k * pow(k / (omega * k * sqrt(1 / r)), abs(l)) * laguerre(p, abs(l), 2 * pow(r, 2) * pow(k, 2) / (pow(omega*k, 2)));
			//cout << phasorOut << endl;
			
			double realField = real(phasorOut);
			double imagField = imag(phasorOut);

			ReEField.at(i).at(j).at(0) = realField;
			ImEField.at(i).at(j).at(0) = imagField;
			
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

double GaussianBeam::laguerre(double k, double alpha, double x) {
	if (k == 0) return 1;
	if (k == 1) return 1 + alpha - x;

	return ((2 * k + 1 + alpha - x)*laguerre(x, alpha, k) - (k + alpha)*laguerre(x, alpha, k - 1)) / (k + 1);
}

int GaussianBeam::getDims() {
	return dimset;
}

double GaussianBeam::realEAt(int i, int j, int k) {
	return ReEField.at(i).at(j).at(k);
}

double GaussianBeam::imagEAt(int i, int j, int k) {
	return ImEField.at(i).at(j).at(k);
}

int GaussianBeam::sign(double val) {
	if (val > 0) return 1;
	if (val < 0) return -1;
	return 0;
}

void GaussianBeam::rootGraph_3d(int argc, char** argv, vector<vector<vector<double> > > Field){
//    //Open root graphics
//    TApplication theApp("App", &argc, argv);
//    gStyle->SetOptStat(0);
////    gStyle->SetPalette(82);
//    TCanvas *c1 = new TCanvas("c1","c1",600,600);
//    TH2F *hcontz = new TH2F("hcontz","Gaussian Beam Cross Section",40,xMin,xMax,40,yMin,yMax);
//    Float_t px, py;
//	for (int i = 0; i < xVals.size(); ++i) {
//		for (int j = 0; j < yVals.size(); ++j) {
//				hcontz->Fill(xVals.at(i),yVals.at(j),Field.at(i).at(j).at(0));
//		}
//	}
////    hcontz->SetTitle("Cross Section of Gaussian Beam");
//    hcontz->GetXaxis()->SetTitle("x"); 
//    hcontz->GetYaxis()->SetTitle("y"); 
//    hcontz->SetMarkerStyle(1);
//    hcontz->GetXaxis()->CenterTitle(); 
//    hcontz->GetYaxis()->CenterTitle(); 
//    hcontz->SetMaximum(4);
//    hcontz->SetMinimum(0);    
//    hcontz->Draw("CONT1z");
//    // Output PDF
//    c1->Print("GBplots.pdf","pdf");
//    theApp.Run();
}

void GaussianBeam::rootGraph_2d(int argc, char** argv, Int_t dim, vector<vector<double > > Points){
    TApplication theApp("App", &argc, argv);
   TCanvas *c1 = new TCanvas("c1","c1",200,10,600,400);
   c1->SetFillColor(42);
   c1->SetGrid();
   const Int_t n = 20;
   Double_t x[n], y[n];
   for (Int_t i=0;i<n;i++) {
      x[i] = i*0.1;
      y[i] = 10*sin(x[i]+0.2);
   }
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerSize(1.5);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Option ACP example");
   gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gr->Draw("ACP");
   // TCanvas::Update() draws the frame, after which one can change it
}

//void GaussianBeam::rootGraph_2d(int argc, char** argv, Int_t dim, vector<vector<double > > Points){
//    TApplication theApp("App", &argc, argv);
////    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
//    TCanvas *c1 = new TCanvas("c1","c1",600,600);
//    Double_t x[100], y[100];
//    Int_t n = dim;
//    for (Int_t i=0;i<n;i++) {
//      x[i] = Points.at(i).at(0);
//      y[i] = Points.at(i).at(1);
//    }
//    TGraph* gr = new TGraph(n,x,y);
//    gr->Draw("AC*");
//}
