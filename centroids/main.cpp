#include "GaussianBeam.h"

const static double PI = 3.14159265;
const static double NGLASS = 1.3;
const static double THETA = 30.0;
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

vector<complex<double> > ETildeBase (vector<complex<double> > f, double theta, vector<complex<double> > kVecs) {
	vector<complex<double> > theVec(3, complex<double> (0,0));

	theVec.at(0) = f.at(0)*refinTM(NGLASS, theta) - f.at(1)*kVecs.at(0)*(1 / tan(theta*PI/180))*(refinTM(NGLASS, theta) + refinTE(NGLASS, theta));
	theVec.at(1) = f.at(1)*refinTE(NGLASS, theta) + f.at(0)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refinTM(NGLASS, theta) + refinTE(NGLASS, theta));
	theVec.at(2) = -f.at(0)*refinTM(NGLASS, theta)*kVecs.at(0) - f.at(1)*refinTE(NGLASS, theta)*kVecs.at(1);

	return theVec;
}

int main(int argc, char** argv){
    GaussianBeam beam1(1.0,.01,1,0);
    beam1.calculateGaussData();
    //If we're only worrying about the interface we need to pick option 2
    //(Basically a 2d, 3d vector)

    //Generate the input and output vectors for fftw
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * beam1.getRealE().size() * beam1.getRealE().at(0).size());
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * beam1.getRealE().size() * beam1.getRealE().at(0).size());

    //Create plane for forward transform
	fftw_plan g = fftw_plan_dft_2d(beam1.getRealE().size(), beam1.getRealE().at(0).size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	int k = 0;
	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			in[k][0] = beam1.getRealE().at(i).at(j).at(0);
			in[k][1] = beam1.getImE().at(i).at(j).at(0);
			k++;
		}
	}

    //Do forward forier transform. Results are stored in "out"
	fftw_execute(g);

    //Seperate the real and imaginart parts of the fourier data (I don't know if
    //we need to do this)	
    vector<vector<vector<double> > > ReFour(beam1.getRealE().size(), vector<vector<double> >(beam1.getRealE().at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<double> > > ImFour(beam1.getRealE().size(), vector<vector<double> >(beam1.getRealE().at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<complex<double> > > > FourData(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(1, complex<double>(0, 0))));

	//k = 0;
	//for (int i = 0; i < beam1.getRealE().size(); i++) {
	//	for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
	//		ReFour.at(i).at(j).at(0) = out[k][0];
	//		ImFour.at(i).at(j).at(0) = out[k][1];
	//		k++;
	//	}
	//}

	k = 0;
	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			complex<double> fourEnt(out[k][0], out[k][1]);
			FourData.at(i).at(j).at(0) = fourEnt;
			k++;
		}
	}

    //Define reflection coefficients for TE and TM
    double refTE = refinTE(NGLASS, THETA);
	double refTM = refinTM(NGLASS, THETA);

	//Assuming horizontal polarization. According to Centroid Shifts paper, f={1,0,0}
    vector<complex<double> > fVec(3, complex<double>(0, 0));
    fVec.at(0)=1;
    fVec.at(1)=0;
    fVec.at(2)=0;

    //Generate the kapa table (using the step size found in the mathematica nb "Single interface Shifts")
	
    //vector<double> kXVals, kYVals; This Declaration wasn't being used... Don't know if we need it.
	vector<vector<vector<complex<double> > > > kPerpTab(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(3, complex<double>(0, 0))));

	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			kPerpTab.at(i).at(j).at(0) = generateK(i, beam1.getRealE().size(), beam1.getK(), findMax(beam1.getXVals()));
			kPerpTab.at(i).at(j).at(1) = generateK(j, beam1.getRealE().at(0).size(), beam1.getK(), findMax(beam1.getYVals()));
			kPerpTab.at(i).at(j).at(2) = sqrt(1-pow(generateK(j, beam1.getRealE().at(0).size(), beam1.getK(), findMax(beam1.getYVals())),2)- pow(generateK(i, beam1.getRealE().size(), beam1.getK(), findMax(beam1.getXVals())), 2));
		}
	}

    //From the Mathematica code:
    //eRtab = Table[eR /. {\[Kappa]x -> \[Kappa]tab[[i, j]][[1]], \[Kappa]y -> \[Kappa]tab[[i, j]][[2]]} /. params$here, {i, 1, dimset}, {j, 1, dimset}];
    //we generate the eR table.
    //NOTE: May be able to compress

	vector<vector<vector<complex<double> > > > eRTab(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(3, complex<double>(0, 0))));

	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			eRTab.at(i).at(j).at(0) = ETildeBase(fVec, THETA, kPerpTab.at(i).at(j)).at(0);
			eRTab.at(i).at(j).at(1) = ETildeBase(fVec, THETA, kPerpTab.at(i).at(j)).at(1);
			eRTab.at(i).at(j).at(2) = ETildeBase(fVec, THETA, kPerpTab.at(i).at(j)).at(2);
		}
	}
	
    //Now we can use the eRtab and the fourier data "out" (E$Tilde in the
    //mathematica notebook) to generate ERtab

	vector<vector<vector<complex<double> > > > ERTab(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(3, complex<double>(0, 0))));

/*	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			ERTab.at(i).at(j).at(0) = eRTab.at(i).at(j).at(0) * //Neet a complex vector that contains fourier data to multiply by 
			ERTab.at(i).at(j).at(1) = eRTab.at(i).at(j).at(1) * //Neet a complex vector that contains fourier data to multiply by
			ERTab.at(i).at(j).at(2) = eRTab.at(i).at(j).at(2) * //Neet a complex vector that contains fourier data to multiply by
		}
	}
 */       
    fftw_plan h = fftw_plan_dft_2d(beam1.getRealE().size(), beam1.getRealE().at(0).size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	
//    beam1.rootGraph(argc, argv, 0);
    return 0;
}
