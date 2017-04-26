#include "GaussianBeam.h"

// Looking for shift of -1.08401 @ 44 degrees and n=0.6592

const static double PI = 3.14159265;
const static double NVAL = 0.659283;
const static double THETA = 60.0;
using namespace std;

complex<double> rTE(double n, double theta, vector<double> kvec){
	theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
	double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);
	//cout << theDot << endl;

    double beta= acos(theDot);
	//cout << beta << endl;
	complex<double> arg(pow(n, 2) - pow(sin(beta), 2), 0);
	complex<double> ans = (cos(beta) - sqrt(arg)) / (cos(beta) + sqrt(arg));
	/*if (isnan(ans.real())) {
		complex<double> backup(0, 0);
		return backup;
	}*/
	return ans;
} 

complex<double> rTM(double n, double theta, vector<double> kvec){
	theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
	double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);
	//cout << theDot << endl;

	double beta = acos(theDot);
	//cout << beta << endl;
	complex<double> arg(pow(n, 2) - pow(sin(beta), 2), 0);
	complex<double> ans = (pow(n,2)*cos(beta)-sqrt(arg))/(pow(n,2)*cos(beta)+sqrt(arg));
	/*if (isnan(ans.real())) {
		complex<double> backup(0, 0);
		return backup;
	}*/
	return ans;
} 

double generateK(int index, int dimsize, int k, int xmax) {
	return (2 * PI) / (2 * k*xmax)*(index-(dimsize+1)/2);
}

double findMax(vector<double> vals) {
	double max = vals.at(0);
	for (int i = 1; i < vals.size(); i++) if (vals.at(i) >= max) max = vals.at(i);
	return max;
}

vector<complex<double> > eRBase (vector<complex<double> > f, double theta, vector<double>  REkVecs) {
	vector<complex<double> > theVec(3, complex<double> (0,0));
	vector<complex<double> > kVecs(0, complex<double> (0,0));
    kVecs.push_back(complex<double>(REkVecs.at(0),0));
    kVecs.push_back(complex<double>(REkVecs.at(1),0));
    kVecs.push_back(complex<double>(REkVecs.at(2),0));
    complex<double> refTM = rTM(NVAL, theta, REkVecs);
    complex<double> refTE = rTE(NVAL, theta, REkVecs);

	theVec.at(0) = f.at(0)*refTM - f.at(1)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
	theVec.at(1) = f.at(1)*refTE + f.at(0)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
	theVec.at(2) = -f.at(0)*refTM*kVecs.at(0) - f.at(1)*refTE*kVecs.at(1);
    /*cout << "(" << real(rTM(NVAL, theta, REkVecs)) << "," << imag(rTM(NVAL, theta, REkVecs)) << ")\t(" << real(rTE(NVAL, theta, REkVecs)) << "," <<  imag(rTE(NVAL, theta, REkVecs)) << ")" << endl;
	char dum;
	cin >> dum;
	*/
	return theVec;
}

double chop(double num) {
	if (abs(num) < 1e-12) return 0;
	return num;
}

int main(int argc, char** argv){
	time_t start = clock();
    cout << "Starting Calculations" << endl;
    double k0 = 1;
    GaussianBeam beam1(20/k0,2*PI,0,0);
    beam1.calculateGaussData();
	int dimset = beam1.getDims();
	//cout << "<" << beam1.getRealE().at(0).at(0).at(0) << " ," << beam1.getImE().at(0).at(0).at(0) << ">" << endl;

	//Assuming horizontal polarization. According to Centroid Shifts paper, f={1,0,0}
    vector<complex<double> > fVec(3, complex<double>(0, 0));
    fVec.at(0)=1;
    fVec.at(1)=0;
    fVec.at(2)=0;

	ofstream fout;
	fout.open("waveTest.txt");
	if (fout.fail()) {
		cerr << "Something went wrong" << endl;
		exit(1);
	}

    //Define the max xKappa and yKappa values
	double xKappa = findMax(beam1.getXVals());
	double yKappa = findMax(beam1.getYVals());	
	
    //Generate the input and output vectors for fftw
	fftw_complex *in, *out, *inx, *outx, *iny, *outy, *inz, *outz, *in3, *out3;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimset * dimset);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dimset * dimset);

    //Create plane for forward transform
	fftw_plan g = fftw_plan_dft_2d(dimset, dimset, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	int k = 0;
	for (int i = 0; i < dimset; i++) {
		for (int j = 0; j < dimset; j++) {
			in[k][0] = beam1.realEAt(i, j, 0)*exp(-100);
			in[k][1] = beam1.imagEAt(i, j, 0)*exp(-100);
            //cout << "<" << chop(beam1.getRealE().at(i).at(j).at(0)) << "," << chop(beam1.getImE().at(i).at(j).at(0)) << ">" << "\t"; 
            //cout << "<" << chop(beam1.getRealE().at(i).at(j).at(0)*exp(-100)) << "," << chop(beam1.getImE().at(i).at(j).at(0)*exp(-100)) << ">" << "\t"; 
            //cout << beam1.getImE().at(i).at(j).at(0) << endl;
			k++;
			//cout << k++ << endl;
		}
       // cout << endl;
	}

	cout << "Checkpoint: Fourier array prepared. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;
    //Do forward forier transform. Results are stored in "out"
	fftw_execute(g);

	vector<vector<vector<complex<double> > > > FourData(dimset, vector<vector<complex<double> > > (dimset, vector<complex<double> >(1, complex<double>(0, 0))));

	k = 0;
	for (int i = 0; i < dimset; i++) {
		for (int j = 0; j < dimset; j++) {
			complex<double> fourEnt(out[k][0] / dimset, out[k][1] / dimset);
			FourData.at(i).at(j).at(0) = fourEnt;
            if (pow(generateK(i, dimset, beam1.getK(), xKappa),2) + pow(generateK(j, dimset, beam1.getK(), yKappa),2) >= 1.0) 
				FourData.at(i).at(j).at(0) = (0,0); //Remove evanescence
            //cout << FourData.at(i).at(j).at(0) << endl;
			k++;
		}
	}

	for (int i = 0; i < FourData.size(); i++) {
		for (int j = 0; j < (FourData.size() + 1) / 2; j++) {
			FourData.at(i).push_back(FourData.at(i).at(0));
			FourData.at(i).erase(FourData.at(i).begin());
		}
	}

	for (int i = 0; i < (FourData.size() + 1) / 2; i++) {
		FourData.push_back(FourData.at(0));
		FourData.erase(FourData.begin());
	}

	//cout << FourData.at(0).at(0).at(0) << endl;
	
	//for (int i = 0; i < FourData.size(); i++) for (int j = 0; j < FourData.at(0).size(); j++) cout << FourData.at(i).at(j).at(0) << endl;

	//cout << "Checkpoint: Fourier data created and rotated. " << (clock() - start) / CLOCKS_PER_SEC << " seconds." << endl;

    //From the Mathematica code:
    //eRtab = Table[eR /. {\[Kappa]x -> \[Kappa]tab[[i, j]][[1]], \[Kappa]y -> \[Kappa]tab[[i, j]][[2]]} /. params$here, {i, 1, dimset}, {j, 1, dimset}];
    //we generate the eR table.
    //NOTE: May be able to compress

	vector<vector<vector<complex<double> > > > eRTab(dimset, vector<vector<complex<double> > >(dimset, vector<complex<double> >(3, complex<double>(0, 0))));

	for (int j = 0; j < dimset; j++) {
		for (int i = 0; i < dimset; i++) {
			vector<double> kVec = vector<double> (0,0);
			double kx = generateK(i + 1, dimset, beam1.getK(), xKappa);
			double ky = generateK(j + 1, dimset, beam1.getK(), yKappa);

			kVec.push_back(kx);
			kVec.push_back(ky);

			double kz = sqrt(1 - pow(kx, 2) - pow(ky, 2));
			if (isnan(kz)) kVec.push_back(0);
			else kVec.push_back(kz);
			eRTab.at(i).at(j) = eRBase(fVec, THETA, kVec);
			/*cout << "<" << kVec.at(0) << "," << kVec.at(1) << "," << kVec.at(2) << ">" << endl;
			char dum;
			cin >> dum;
			*/
			//cout << "<" << eRTab.at(i).at(j).at(0).real() << "," << eRTab.at(i).at(j).at(0).imag() << ">" << endl;
		}
	}

	//cout << "<" << eRTab.at(0).at(0).at(0) << "," << eRTab.at(0).at(0).at(1) << "," << eRTab.at(0).at(0).at(2) << ">" << endl;
	//for (int j = 0; j < beam1.getRealE().size(); j++) for (int i = 0; i < beam1.getRealE().at(0).size(); i++)	cout << "<" << eRTab.at(i).at(j).at(0) << "," << eRTab.at(i).at(j).at(1) << "," << eRTab.at(i).at(j).at(2) << ">" << endl;
	
	cout << "Checkpoint: ETilde generated. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;
    //Now we can use the eRtab and the fourier data "out" (E$Tilde in the
    //mathematica notebook) to generate ERtab

	vector<vector<vector<complex<double> > > > ERTab(dimset, vector<vector<complex<double> > >(dimset, vector<complex<double> >(3, complex<double>(0, 0))));
    int nanc = 0;

    for (int i = 0; i < dimset; i++) {
    	for (int j = 0; j < dimset; j++) {
			complex<double> fourPoint(FourData.at(i).at(j).at(0).real(), -FourData.at(i).at(j).at(0).imag());

    		ERTab.at(i).at(j).at(0) = eRTab.at(i).at(j).at(0) * fourPoint;
    		ERTab.at(i).at(j).at(1) = eRTab.at(i).at(j).at(1) * fourPoint;
    		ERTab.at(i).at(j).at(2) = eRTab.at(i).at(j).at(2) * fourPoint;

			if (isnan(ERTab.at(i).at(j).at(0).real()) || isnan(ERTab.at(i).at(j).at(1).real()) || isnan(ERTab.at(i).at(j).at(2).real()) || isnan(ERTab.at(i).at(j).at(0).imag()) || isnan(ERTab.at(i).at(j).at(1).imag()) || isnan(ERTab.at(i).at(j).at(2).imag()))
            cout << ++nanc << " NAN found: " << i << ", " << j << endl;
                
            //cout << ++nanc << " NAN found: " << i << ", " << j << endl;
            
			/*cout << "<(" << chop(real(ERTab.at(i).at(j).at(0))) << "," << chop(imag(ERTab.at(i).at(j).at(0))) << "), (" << chop(real(ERTab.at(i).at(j).at(1))) << "," << chop(imag(ERTab.at(i).at(j).at(1))) << "), (" << chop(real(ERTab.at(i).at(j).at(2))) << "," << chop(imag(ERTab.at(i).at(j).at(2))) << ")>" << endl;
			char dum;
			cin >> dum;
			*/
			}
    }

	//for (int i = 0; i < 30; i++) cout << "<(" << real(ERTab.at(0).at(i).at(0)) << "," << imag(ERTab.at(0).at(i).at(0)) << "), (" << real(ERTab.at(0).at(i).at(1)) << "," << imag(ERTab.at(0).at(i).at(1)) << "), (" << real(ERTab.at(0).at(i).at(2)) << "," << imag(ERTab.at(0).at(i).at(2)) << ")>" << endl;

    //cout << "<(" << chop(real(ERTab.at(5).at(5).at(0))) << "," << chop(imag(ERTab.at(5).at(5).at(0))) << "), (" << chop(real(ERTab.at(5).at(5).at(1))) << "," << chop(imag(ERTab.at(5).at(5).at(1))) << "), (" << chop(real(ERTab.at(5).at(5).at(2))) << "," << chop(imag(ERTab.at(5).at(5).at(2))) << ")>" << endl;
	
	//DO NOT CHANGE ANYTHING ABOVE THIS LINE!! AS FAR AS WE KNOW, IT'S WORKING FINE!
	cout << "Checkpoint: Ready for IFFT. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

	inx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
	outx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
    
	iny = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
	outy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());

	inz = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());
	outz = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size());

	in3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);
	out3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);

    k = 0;
    for (int i = 0; i < ERTab.size(); i++) {
        for (int j = 0; j < ERTab.at(0).size(); j++){ 
			inx[k][0] = real(ERTab.at(i).at(j).at(0));
			inx[k][1] = -imag(ERTab.at(i).at(j).at(0));

			iny[k][0] = real(ERTab.at(i).at(j).at(1));
			iny[k][1] = -imag(ERTab.at(i).at(j).at(1));

			inz[k][0] = real(ERTab.at(i).at(j).at(2));
			inz[k][1] = -imag(ERTab.at(i).at(j).at(2));
			k++;
        }    
    }

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
			in3[k][0] = real(ERTab.at(i).at(j).at(0));
			in3[k][1] = -imag(ERTab.at(i).at(j).at(0));
			k++;

			in3[k][0] = real(ERTab.at(i).at(j).at(1));
			in3[k][1] = -imag(ERTab.at(i).at(j).at(1));
			k++;

			in3[k][0] = real(ERTab.at(i).at(j).at(2));
			in3[k][1] = -imag(ERTab.at(i).at(j).at(2));
			k++;
		}
	}

    fftw_plan hx = fftw_plan_dft_2d(ERTab.size(), ERTab.at(0).size(), inx, outx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(hx);

	fftw_plan hy = fftw_plan_dft_2d(ERTab.size(), ERTab.at(0).size(), iny, outy, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(hy);

	fftw_plan hz = fftw_plan_dft_2d(ERTab.size(), ERTab.at(0).size(), inz, outz, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(hz);

	fftw_plan h3 = fftw_plan_dft_3d(ERTab.size(), ERTab.at(0).size(), 3, in3, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(h3);

	vector<vector<vector<complex<double> > > > outBeam(ERTab.size(), vector<vector<complex<double> > >(ERTab.at(0).size(), vector<complex<double> > (3, complex<double> (0,0))));
    vector<vector<vector<double> > > REoutBeam(ERTab.size(), vector<vector<double> > (ERTab.at(0).size(), vector<double> (3,0)));
	vector<vector<vector<double> > > outBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<double> > > OGBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));

	//cout << "After IFT:" << endl;

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
				vector<complex<double> > fourEnt;
				/*complex<double> entry1(outx[k][0] / ERTab.size(), outx[k][1] / ERTab.size());
				fourEnt.push_back(entry1);

				complex<double> entry2(outy[k][0] / ERTab.size(), outy[k][1] / ERTab.size());
				fourEnt.push_back(entry2);

				complex<double> entry3(outz[k][0] / ERTab.size(), outz[k][1] / ERTab.size());
				fourEnt.push_back(entry3);
				k++;
				*/
				complex<double> entry1(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
				fourEnt.push_back(entry1);
				k++;

				complex<double> entry2(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
				fourEnt.push_back(entry2);
				k++;

				complex<double> entry3(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
				fourEnt.push_back(entry3);
				k++;

				outBeam.at(i).at(j) = fourEnt;
				//cout << "<" << fourEnt.at(0) << ", " << fourEnt.at(1) << ", " << fourEnt.at(2) << ">" << endl;
				//cout << "<" << chop(fourEnt.at(0).real()) << ", " << chop(fourEnt.at(0).imag()) << ", " << ">" << "\t";
				//cout << "<" << chop(sqrt(pow(fourEnt.at(0).real(),2) + pow(fourEnt.at(0).imag(),2))) << ">" << "\t";
		}
        //cout << endl;
	}

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
			REoutBeam.at(i).at(j).at(0) = outx[k][0];
			REoutBeam.at(i).at(j).at(1) = outy[k][0];
			REoutBeam.at(i).at(j).at(2) = outz[k][0];
			k++;
		}
	}

	cout << "Checkpoint: Heavy lifting complete. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

    fftw_free(in);
    fftw_free(out);
    fftw_free(inx);
    fftw_free(outx);
	fftw_free(iny);
	fftw_free(outy);
	fftw_free(inz);
	fftw_free(outz);
	fftw_free(in3);
	fftw_free(out3);

    fftw_destroy_plan(hx);
	fftw_destroy_plan(hy);
	fftw_destroy_plan(hz);
	fftw_destroy_plan(h3);
    fftw_destroy_plan(g);

	fout.close();

    //Calculate magnitude of both beams for comparision
	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		outBeamMag.at(i).at(j).at(0) = sqrt(real(outBeam.at(i).at(j).at(0) * conj(outBeam.at(i).at(j).at(0))) + real(outBeam.at(i).at(j).at(1) * conj(outBeam.at(i).at(j).at(1))) + real(outBeam.at(i).at(j).at(2) * conj(outBeam.at(i).at(j).at(2))));

	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		OGBeamMag.at(i).at(j).at(0) = sqrt(pow(beam1.realEAt(i,j,0),2) + pow(beam1.imagEAt(i, j, 0),2));

    //Calculate centroid Shifts
	double nXrp1 = 0, nYrp1 = 0, denom = 0;
	
	for (int i = 0; i < outBeamMag.size(); i++) {
		for (int j = 0; j < outBeamMag.at(0).size(); j++) {
				//cout << (i + 1) * pow(outBeamMag.at(j).at(i).at(0), 2) << endl;
				nXrp1 += (j + 1) * pow(outBeamMag.at(j).at(i).at(0), 2);
				nYrp1 += (i + 1) * pow(outBeamMag.at(j).at(i).at(0), 2);
				denom += pow(outBeamMag.at(i).at(j).at(0), 2);
				//cout << nXrp1 << " " << nYrp1 << " " << denom << endl;
		}
	}

	double nXr = nXrp1 / denom;
	double nYr = nYrp1 / denom;

    //Compare calculated shifts to analytical result
    complex<double> ARshift1 = (4*pow(NVAL,2)*sin(THETA*PI/180))/(beam1.getK()*(-1+pow(NVAL,2)+(1+pow(NVAL,2))*cos(2*THETA*PI/180))*sqrt(complex<double>(-1*pow(NVAL,2)+pow(sin(THETA*PI/180),2),0)));
    complex<double> ARshift2 = -2*sqrt(2)*sin(THETA*PI/180)/(beam1.getK()*sqrt(complex<double>(1-2*pow(NVAL,2)-cos(2*THETA*PI/180),0)));

    cout << "Calculated" << endl;
    cout << "(" << nXr << "," << nYr << ")" << endl;

    //cout << "Analytical" << endl;
    //cout << "(" << ARshift1 << "," << ARshift2 << ")" << endl;

	cout << "Total time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

   // beam1.rootGraph(argc, argv, OGBeamMag);
    //beam1.rootGraph(argc, argv, outBeamMag);
    return 0;
}
