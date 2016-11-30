#include "GaussianBeam.h"

#include<complex>
#include<cmath>

const static double PI = 3.14159265;
const static double NGLASS = 1.3;
using namespace std;

double snell(double n, double thetaIrad) {
	double arg, thetaTrad;
	arg = 1 / n * sin(thetaIrad);
	thetaTrad = asin(arg);
	return thetaTrad;
}

double refcofTE(double n, double thetaI) {
	double thetaTrad, thetaIrad, refcofTE;
	thetaIrad = thetaI*PI / 180;
	thetaTrad = snell(n, thetaIrad);
	refcofTE = abs((cos(thetaIrad) - n*cos(thetaTrad)) / (cos(thetaIrad) + n*cos(thetaTrad)));
	return refcofTE;
}

double refinTE(double n, double thetaI) {
	double thetaTrad, thetaIrad, refcofTE, refinTE;
	thetaIrad = thetaI*PI / 180;
	thetaTrad = snell(n, thetaIrad);
	refcofTE = abs((cos(thetaIrad) - n*cos(thetaTrad)) / (cos(thetaIrad) + n*cos(thetaTrad)));
	refinTE = pow(refcofTE, 2);
	return refinTE;
}

double refcofTM(double n, double thetaI) {
	double thetaTrad, thetaIrad, refcofTM;
	thetaIrad = thetaI*PI / 180;
	thetaTrad = snell(n, thetaIrad);
	refcofTM = abs((cos(thetaTrad) - n*cos(thetaIrad)) / (cos(thetaTrad) + n*cos(thetaIrad)));
	return refcofTM;
}

double refinTM(double n, double thetaI) {
	double thetaTrad, thetaIrad, refcofTM, refinTM;
	thetaIrad = thetaI*PI / 180;
	thetaTrad = snell(n, thetaIrad);
	refcofTM = abs((cos(thetaTrad) - n*cos(thetaIrad)) / (cos(thetaTrad) + n*cos(thetaIrad)));
	refinTM = pow(refcofTM, 2);
	return refinTM;
}

double generateK(int index, int dimsize, int k, int xmax) {
	return (2 * PI) / (2 * k*xmax)*(index-(dimsize+1)/2);
}

double findMax(vector<double> vals) {
	double max = vals.at(0);
	for (int i = 1; i < vals.size(); i++) if (vals.at(i) >= max) max = vals.at(i);
	return max;
}

vector<complex<double> > ETildeBase (vector<complex<int> > f, double theta, vector<complex<double> > kVecs) {
	vector<complex<double> > theVec(3, complex<double> (0,0));

	theVec.at(0) = f.at(0)*refinTM(NGLASS, theta) - f.at(1)*kVecs.at(0)*cot(theta)*(refinTM(NGLASS, theta) + refinTE(NGLASS, theta));
	theVec.at(1) = f.at(1)*refinTE(NGLASS, theta) + f.at(0)*kVecs.at(1)*cot(theta)*(refinTM(NGLASS, theta) + refinTE(NGLASS, theta));
	theVec.at(2) = -f.at(0)*refinTM(NGLASS, theta)*kVecs.at(0) - f.at(1)*refinTE(NGLASS, theta)*kVecs.at(1);

	return theVec;
}

int main(int argc, char** argv){
    GaussianBeam beam1(1.0,.01,1,0);
    beam1.calculateGaussData();
    //If we're only worrying about the interface we need to pick option 2
    //Flatten to 2D at a specific z   

	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * beam1.getRealE().size() * beam1.getRealE().at(0).size());
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * beam1.getRealE().size() * beam1.getRealE().at(0).size());

	fftw_plan g = fftw_plan_dft_2d(beam1.getRealE().size(), beam1.getRealE().at(0).size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	int k = 0;
	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			in[k][0] = beam1.getRealE().at(i).at(j).at(0);
			in[k][1] = beam1.getImE().at(i).at(j).at(0);
			k++;
		}
	}

	fftw_execute(g);

	vector<vector<vector<double> > > ReFour(beam1.getRealE.size(), vector<vector<double> >(beam1.getRealE.at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<double> > > ImFour(beam1.getRealE.size(), vector<vector<double> >(beam1.getRealE.at(0).size(), vector<double>(1, 0)));

	k = 0;
	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			ReFour.at(i).at(j).at(0) = out[k][0];
			ImFour.at(i).at(j).at(0) = out[k][1];
			k++;
		}
	}

	double refTE = refinTE(1.3, 30);
	double refTM = refinTM(1.3, 30);

	//Assuming horizontal polarization. According to Centroid Shifts paper, f1 = 1, f2 = 0

	vector<double> kXVals, kYVals;
	vector<vector<vector<complex<double> > > > kPerpTab(beam1.getRealE.size(), vector<vector<complex<double> > >(beam1.getRealE.at(0).size(), vector<complex<double> >(3, complex<double>(0, 0)));

	for (int i = 0; i < beam1.getRealE.size(); i++) {
		for (int j = 0; j < beam1.getRealE.at(0).size(); j++) {
			kPerpTab.at(i).at(j).at(0) = generateK(i, beam1.getRealE.size(), beam1.getK(), findMax(beam1.getXVals()));
			kPerpTab.at(i).at(j).at(1) = generateK(j, beam1.getRealE.at(0).size(), beam1.getK(), findMax(beam1.getYVals()));
			kPerpTab.at(i).at(j).at(2) = sqrt(1-pow(generateK(j, beam1.getRealE.at(0).size(), beam1.getK(), findMax(beam1.getYVals())),2)- pow(generateK(i, beam1.getRealE.size(), beam1.getK(), findMax(beam1.getXVals())), 2));
		}
	}



	fftw_plan h = fftw_plan_dft_2d(beam1.getRealE().size(), beam1.getRealE().at(0).size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	
    beam1.rootGraph(argc, argv, 0);
    return 0;
}
