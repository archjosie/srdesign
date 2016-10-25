#include "GaussianBeam.h"
const static double PI = 3.14159265;
using namespace std;

GaussianBeam::GaussianBeam() {
	lambda = 1;
	w0 = 20;
	m = 0;
	n = 0;

	k = 2 * PI / lambda;
	zR = calculateRayleigh();
}

GaussianBeam::GaussianBeam(double w0, double lambda, unsigned int m, unsigned int n) {
	this->w0 = w0;
	this->lambda = lambda;
	this->m = m;
	this->n = n;

	this->k = 2 * PI / lambda;
	this->zR = calculateRayleigh();
}

double GaussianBeam::distance(double x, double y) {
	return sqrt(pow(x, 2) + pow(y, 2));
}

double GaussianBeam::calculateRadCurv(double z) {
	return (z*(1 + pow(zR / z, 2)));
}

double GaussianBeam::calculateGouy(double z) {
	return atan(z / zR);
}

double GaussianBeam::calculateWaist(double z) {
	return w0*sqrt(1 + pow(z / zR, 2));
}

double GaussianBeam::calculateRayleigh() {
	return PI * pow(w0, 2) / lambda;
}

double GaussianBeam::calculateHermite(double x, unsigned int m) { //Calculates the nth Hermite polynomial evaluated at x recursively
	if (m == 0) return 1;
	if (m == 1) return 2 * x; //Base cases

	return (2 * x*calculateHermite(x, m - 1) - 2 * (m - 1)*calculateHermite(x, m - 2));
}

void GaussianBeam::calculateGaussData() {
	bool choiceMade = false;
	int choice = 0;

	do {
		cout << "Please tell us what type of data you would like us to make: " << endl;
		cout << "1: Gaussian beam data across range of x,y,z" << endl;
		cout << "2: Gaussian beam data at fixed z (cross-sectional snapshot)" << endl;
		cout << "3: Fourier transform data of Gaussian beam across cross-section." << endl;

		cin >> choice;

		if (choice > 0 && choice < 4) {
			choiceMade = true;
		}
		else {
			cout << "Invalid selection. Please try again." << endl << endl;
		}
	} while (!choiceMade);

	double xMax; //Stores upper limit on x
	cout << "Enter the maximum value of x: " << endl;
	cin >> xMax;

	double xMin; //Stores upper limit on radius of consideration
	cout << "Enter the minimum value of x: " << endl;
	cin >> xMin;

	double xInt; //Stores number of radial steps in the data output
	cout << "Enter the step size on x: " << endl;
	cin >> xInt;

	double yMax; //Stores upper limit on radius of consideration
	cout << "Enter the maximum value of y: " << endl;
	cin >> yMax;

	double yMin; //Stores upper limit on radius of consideration
	cout << "Enter the minimum value of y: " << endl;
	cin >> yMin;

	double yInt; //Stores number of radial steps in the data output
	cout << "Enter the step size on y: " << endl;
	cin >> yInt;

	vector<double> xVals;
	vector<double> yVals;

	int xRange = (xMax - xMin) / xInt + 1;
	int yRange = (yMax - yMin) / yInt + 1;

	for (int i = 0; i < xRange; i++) { //Populate xVals
		double xCurr = xMin + i*xInt;
		xVals.push_back(xCurr);
	}

	for (int i = 0; i < yRange; i++) { //Populate yVals
		double yCurr = yMin + i*yInt;
		yVals.push_back(yCurr);
	}

	vector<double> zVals;

	if (choice == 1) {

		double zMin; //Stores lower limit on distance from focus
		cout << "Enter the minimum distance from focus: " << endl;
		cin >> zMin;

		double zMax; //Stores upper limit on distance from focus
		cout << "Enter the maximum distance from focus: " << endl;
		cin >> zMax;

		double zInt; //Stores number of distance steps in the data output
		cout << "Enter the step size on z: " << endl;
		cin >> zInt;

		int zRange = (zMax - zMin) / zInt + 1; //"Hackfix" to get as close as possible to max without going over


		for (int i = 0; i < zRange; i++) { //Populate zVals
			double zCurr = i*zInt + zMin; //We don't want 0 in zVals
			zVals.push_back(zCurr);
		}	

	}

	else {
		double theZ; //Stores lower limit on distance from focus
		cout << "Enter the distance from focus you'd like to examine: " << endl;
		cin >> theZ;

		zVals.push_back(theZ);
	}

	vector<vector<vector<double> > > ReEField(xVals.size(), vector<vector<double> >(yVals.size(), vector<double>(zVals.size(), 0)));
	vector<vector<vector<double> > > ImEField(xVals.size(), vector<vector<double> >(yVals.size(), vector<double>(zVals.size(), 0))); //Separate 3D vectors to hold real and imaginary parts of E-field

	for (int l = 0; l < zVals.size(); ++l) {  //Weird choice, but this streamlines the program
										  //These only vary with z, so we only need to calculate them once per loop
		double gouy = calculateGouy(zVals.at(l));
		double radCurv = calculateRadCurv(zVals.at(l));
		double spotSize = calculateWaist(zVals.at(l));

		for (int i = 0; i < xVals.size(); ++i) {
			for (int j = 0; j < yVals.size(); ++j) {

				double r = distance(xVals.at(i), yVals.at(j));
				//Complex argument for Gaussian E-field

				double imArg = -(k* zVals.at(l) + k*(pow(r, 2) / (radCurv * 2)) - gouy);
				complex<double> expArg(0.0, imArg);

				complex<double> phasorOut; //Uncomment if we only want to consider real component of field			
				phasorOut = calculateHermite(sqrt(2)*xVals.at(i) / spotSize, m)*calculateHermite(sqrt(2)*yVals.at(j) / spotSize, n)*(w0 / spotSize)*exp(-pow(r / spotSize, 2))*exp(expArg);
				double realField = real(phasorOut);
				double imagField = imag(phasorOut);

				ReEField.at(i).at(j).at(l) = realField;
				ImEField.at(i).at(j).at(l) = imagField;
			}
		}
	}
	ofstream fout;

	string fileName; //Stores number of distance steps in the data output
	cout << "Data has been generated. Please enter name and format of file to be printed (.dat preferred): " << endl;
	cin >> fileName;

	fout.open(fileName.c_str());

	if (fout.fail()) {
		cerr << "Something went wrong! Exiting now..." << endl;
		exit(1);
	}


	if (choice == 3) {
		fourierTran(ReEField, ImEField, xVals, yVals, fout);
		return;
	}

	//Time to write and export the data
	fout << "# x " << '\t' << "y " << '\t' << "z " << '\t' << "Re(E(x,y,z)) " << '\t' << "Im(E(x,y,z)) " << '\t' << "Abs(E(x,y,z))" << endl; //Pound sign forces gnuplot to ignore line

																																			 //Write our 3D complex vector to file
	for (int i = 0; i < xVals.size(); ++i) {
		for (int j = 0; j < yVals.size(); ++j) {
			for (int k = 0; k < zVals.size(); ++k) {

				//*inhales*
				fout << xVals.at(i) << '\t' << yVals.at(j) << '\t' << zVals.at(k) << '\t' << ReEField.at(i).at(j).at(k) << '\t' << ImEField.at(i).at(j).at(k) << '\t' << distance(ReEField.at(i).at(j).at(k), ImEField.at(i).at(j).at(k)) << endl;
			}
		}
	}

	cout << "Data output complete! Thank you for your business!" << endl;
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
