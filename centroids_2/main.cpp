#include "GaussianBeam.h"

// Looking for shift of -1.08401 @ 44 degrees and n=0.6592

const static double PI = 3.14159265;
const static double NVAL = 0.659283;
const static double THETA = 44.0;
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

vector<complex<double> > ETildeBase (vector<complex<double> > f, double theta, vector<double>  REkVecs) {
	vector<complex<double> > theVec(3, complex<double> (0,0));
	vector<complex<double> > kVecs(0, complex<double> (0,0));
    kVecs.push_back(complex<double>(REkVecs.at(0),0));
    kVecs.push_back(complex<double>(REkVecs.at(1),0));
    kVecs.push_back(complex<double>(REkVecs.at(2),0));

	theVec.at(0) = f.at(0)*refinTM(NVAL, theta) - f.at(1)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refinTM(NVAL, theta) + refinTE(NVAL, theta));
	theVec.at(1) = f.at(1)*refinTE(NVAL, theta) + f.at(0)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refinTM(NVAL, theta) + refinTE(NVAL, theta));
	theVec.at(2) = -f.at(0)*refinTM(NVAL, theta)*kVecs.at(0) - f.at(1)*refinTE(NVAL, theta)*kVecs.at(1);
    cout << refinTM(NVAL, theta) << endl;

	return theVec;
}

int main(int argc, char** argv){
    double k0 = 2*PI/(632.8e-9);
    GaussianBeam beam1(20000/k0,632.8e-9,0,0);
    beam1.calculateGaussData();

    //Define reflection coefficients for TE and TM
    double refTE = refinTE(NVAL, THETA);
	double refTM = refinTM(NVAL, THETA);

	//Assuming horizontal polarization. According to Centroid Shifts paper, f={1,0,0}
    vector<complex<double> > fVec(3, complex<double>(0, 0));
    fVec.at(0)=1;
    fVec.at(1)=0;
    fVec.at(2)=0;

    //Generate the kappa table (using the step size found in the mathematica nb "Single interface Shifts")
	vector<vector<vector<double> > > kPerpTab(beam1.getRealE().size(), vector<vector<double> >(beam1.getRealE().at(0).size(), vector<double> (3, 0)));
	vector<vector<double> > kComp(beam1.getRealE().size(), vector<double>(beam1.getRealE().at(0).size()));

    cout << "We're good up to here" << endl;

	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			//kPerpTab.at(i).at(j).at(0) = generateK(i, beam1.getRealE().size(), beam1.getK(), findMax(beam1.getXVals()));
			//kPerpTab.at(i).at(j).at(1) = generateK(j, beam1.getRealE().at(0).size(), beam1.getK(), findMax(beam1.getYVals()));
			kPerpTab.at(i).at(j).at(2) = sqrt((1-pow(generateK(j, beam1.getRealE().at(0).size(), beam1.getK(), findMax(beam1.getYVals())),2)- pow(generateK(i, beam1.getRealE().size(), beam1.getK(), findMax(beam1.getXVals())), 2)));
			//kComp.at(i).at(j)= pow(generateK(i, beam1.getRealE().size(), beam1.getK(), findMax(beam1.getXVals())),2)+pow(generateK(j, beam1.getRealE().at(0).size(), beam1.getK(), findMax(beam1.getYVals())),2);
		}
	}

    cout << "kappa table generated" << endl;

    //Generate the input and output vectors for fftw
	fftw_complex *in, *out, *in2, *out2;
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

	vector<vector<vector<complex<double> > > > FourData(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(1, complex<double>(0, 0))));

	k = 0;
	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			complex<double> fourEnt(out[k][0], out[k][1]);
			FourData.at(i).at(j).at(0) = fourEnt;
//            if (kComp.at(i).at(j)>=1.0){
//                FourData.at(i).at(j).at(0) = (0,0);
//            }
			k++;
		}
	}

    cout << "Fourier Data Generated" << endl;

    //From the Mathematica code:
    //eRtab = Table[eR /. {\[Kappa]x -> \[Kappa]tab[[i, j]][[1]], \[Kappa]y -> \[Kappa]tab[[i, j]][[2]]} /. params$here, {i, 1, dimset}, {j, 1, dimset}];
    //we generate the eR table.
    //NOTE: May be able to compress

	vector<vector<vector<complex<double> > > > eRTab(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(3, complex<double>(0, 0))));

	for (int i = 0; i < beam1.getRealE().size(); i++) {
		for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
			eRTab.at(i).at(j) = ETildeBase(fVec, THETA, kPerpTab.at(i).at(j));
//            cout << "<" << eRTab.at(i).at(j).at(0)<< "," << eRTab.at(i).at(j).at(1) << "," << eRTab.at(i).at(j).at(2) << ">" << endl;
		}
	}
	
    //Now we can use the eRtab and the fourier data "out" (E$Tilde in the
    //mathematica notebook) to generate ERtab

	vector<vector<vector<complex<double> > > > ERTab(beam1.getRealE().size(), vector<vector<complex<double> > >(beam1.getRealE().at(0).size(), vector<complex<double> >(3, complex<double>(0, 0))));

    for (int i = 0; i < beam1.getRealE().size(); i++) {
    	for (int j = 0; j < beam1.getRealE().at(0).size(); j++) {
    		ERTab.at(i).at(j).at(0) = eRTab.at(i).at(j).at(0) * FourData.at(i).at(j).at(0);
    		ERTab.at(i).at(j).at(1) = eRTab.at(i).at(j).at(1) * FourData.at(i).at(j).at(0);
    		ERTab.at(i).at(j).at(2) = eRTab.at(i).at(j).at(2) * FourData.at(i).at(j).at(0);
    	}
    }
          
    cout << "ERTab generated, process complete" << endl;

	in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
	out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
    
    k = 0;
    for (int i = 0; i < ERTab.size(); i++) {
        for (int j = 0; j < ERTab.at(0).size(); j++) {
            in2[k][0] = ERTab.at(i).at(j).at(0).real();
            in2[k][1] = ERTab.at(i).at(j).at(0).imag();
        }    
    }

    fftw_plan h = fftw_plan_dft_2d(ERTab.size(), ERTab.at(0).size(), in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(h);

	vector<vector<vector<complex<double> > > > outBeam(ERTab.size(), vector<vector<complex<double> > >(ERTab.at(0).size(), vector<complex<double> > (1, complex<double> (0,0))));
    vector<vector<vector<double> > > REoutBeam(ERTab.size(), vector<vector<double> > (ERTab.at(0).size(), vector<double> (1,0)));
	vector<vector<vector<double> > > outBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<double> > > OGBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
			complex<double> fourEnt(out[k][0], out[k][1]);
			outBeam.at(i).at(j).at(0) = fourEnt; 
			k++;
		}
	}

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
			double REfourEnt = out[k][0];
			REoutBeam.at(i).at(j).at(0) = REfourEnt; 
			k++;
		}
	}

    cout << "Output beam constructed." << endl;

    fftw_free(in);
    fftw_free(out);
    fftw_free(in2);
    fftw_free(out2);

    fftw_destroy_plan(h);
    fftw_destroy_plan(g);

    //Calculate magnitude of both beams for comparision
	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		outBeamMag.at(i).at(j).at(0) = sqrt(real(outBeam.at(i).at(j).at(0) * conj(outBeam.at(i).at(j).at(0))));

	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		OGBeamMag.at(i).at(j).at(0) = sqrt(pow(beam1.getRealE().at(i).at(j).at(0),2) + pow(beam1.getImE().at(i).at(j).at(0),2));

    //Calculate centroid Shifts
	double nXrp1 = 0, nYrp1 = 0, denom = 0;
	
	for (int i = 0; i < outBeamMag.size(); i++) {
		for (int j = 0; j < outBeamMag.at(0).size(); j++) {
			nXrp1 += j * outBeamMag.at(i).at(j).at(0);
			nYrp1 += i * outBeamMag.at(i).at(j).at(0);
			denom += outBeamMag.at(i).at(j).at(0);
		}
	}

	double nXr = nXrp1 / denom;
	double nYr = nYrp1 / denom;

    //Compare calculated shifts to analytical result
    complex<double> ARshift1 = (4*pow(NVAL,2)*sin(THETA*PI/180))/(beam1.getK()*(-1+pow(NVAL,2)+(1+pow(NVAL,2))*cos(2*THETA*PI/180))*sqrt(complex<double>(-1*pow(NVAL,2)+pow(sin(THETA*PI/180),2),0)));
    complex<double> ARshift2 = -2*sqrt(2)*sin(THETA*PI/180)/(beam1.getK()*sqrt(complex<double>(1-2*pow(NVAL,2)-cos(2*THETA*PI/180),0)));

    cout << "Calculated" << endl;
    cout << "(" << nXr << "," << nYr << ")" << endl;

    cout << "Analytical" << endl;
    cout << "(" << ARshift1 << "," << ARshift2 << ")" << endl;

//    beam1.rootGraph(argc, argv, OGBeamMag);
    beam1.rootGraph(argc, argv, outBeamMag);
    return 0;
}
