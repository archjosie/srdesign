#include "GaussianBeam.h"

const static double PI = 3.14159265;
const static double NVAL = 0.659283;
using namespace std;

complex<double> rTE(double n, double theta, vector<double> kvec){ //Fresnel coeff for TE waves
//Parameters: refraction index n, incident angle theta (in degrees), E-field vector kvec in Fourier space
	theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
	//Argument for beta: dot product between nvec and kvec
	double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);

    double beta= acos(theDot);
	complex<double> arg(pow(n, 2) - pow(sin(beta), 2), 0);
	complex<double> ans = (cos(beta) - sqrt(arg)) / (cos(beta) + sqrt(arg));
	return ans;
} 

complex<double> rTM(double n, double theta, vector<double> kvec){ //Fresnel coeff for TE waves
//Parameters: refraction index n, incident angle theta (in degrees), E-field vector kvec in Fourier space
	theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
	double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);

	double beta = acos(theDot);
	complex<double> arg(pow(n, 2) - pow(sin(beta), 2), 0);
	complex<double> ans = (pow(n,2)*cos(beta)-sqrt(arg))/(pow(n,2)*cos(beta)+sqrt(arg));
	return ans;
} 

double generateK(int index, int dimsize, int k, int xmax) { //O(1) wave vector generator
	return (2 * PI) / (2 * k*xmax)*(index-(dimsize+1)/2);
}

double findMax(vector<double> vals) { //O(n) maximum search function
	double max = vals.at(0);
	for (int i = 1; i < vals.size(); i++) if (vals.at(i) >= max) max = vals.at(i);
	return max;
}

vector<complex<double> > eRBase (vector<complex<double> > f, double theta, vector<double>  REkVecs) {
	//Models behavior of beam inside interface
	//Parmeters: polarization vector f, incident angle theta (degrees), real-valued wave vector REkVecs
	vector<complex<double> > theVec(3, complex<double> (0,0));
	vector<complex<double> > kVecs(0, complex<double> (0,0));
	
	//Convert RekVecs to complex values
    kVecs.push_back(complex<double>(REkVecs.at(0),0));
    kVecs.push_back(complex<double>(REkVecs.at(1),0));
    kVecs.push_back(complex<double>(REkVecs.at(2),0));
	
	//Calculate Fresnel coefficients
    complex<double> refTM = rTM(NVAL, theta, REkVecs);
    complex<double> refTE = rTE(NVAL, theta, REkVecs);
	
	//Use known beam behavior to return output beam's wave vector
	theVec.at(0) = f.at(0)*refTM - f.at(1)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
	theVec.at(1) = f.at(1)*refTE + f.at(0)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
	theVec.at(2) = -f.at(0)*refTM*kVecs.at(0) - f.at(1)*refTE*kVecs.at(1);
	return theVec;
}

double chop(double num) { //Sets low magnitude numbers to 0
	if (abs(num) < 1e-12) return 0;
	return num;
}

double shift(double THETA){
	time_t start = clock();
    //cout << "Starting Calculations" << endl; //Uncomment for timing studies
    double k0 = 1;
    GaussianBeam beam1(20000/k0,2*PI,0,0); 
    beam1.calculateGaussData(); //Initialize E-field for beam
	int dimset = beam1.getDims();

    //Input three-component polarization vector (no need to normalize). Valid entries are:
	//Horizontal: {1,0,0}
	//Vertical: {0,1,0}
	//CCW Cicular: {1,i,0}
	//CW Circular: {1,-i,0}
    vector<complex<double> > fVec(3, complex<double>(0, 0));
    fVec.at(0)=1;
    fVec.at(1)=complex<double>(0,-1);
    fVec.at(2)=0;

    //Define the max xKappa and yKappa values
	double xKappa = findMax(beam1.getXVals());
	double yKappa = findMax(beam1.getYVals());	
	
    //Allocate memory for input and output arrays for FFTW
	fftw_complex *in, *out, *in3, *out3;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimset * dimset);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dimset * dimset);

    //Create plan for forward transform
	fftw_plan g = fftw_plan_dft_2d(dimset, dimset, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//Load FFTW matrices
	int k = 0;
	for (int i = 0; i < dimset; i++) {
		for (int j = 0; j < dimset; j++) {
			in[k][0] = beam1.realEAt(i, j, 0);
			in[k][1] = beam1.imagEAt(i, j, 0);
			k++;
		}
	}

	//cout << "Checkpoint: Fourier array prepared. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl; //Uncomment for timing studies
	
    //Do forward Fourier transform. Results are stored in "out"
	fftw_execute(g);

	//Transfer results to 3D vector
	vector<vector<vector<complex<double> > > > FourData(dimset, vector<vector<complex<double> > > (dimset, vector<complex<double> >(1, complex<double>(0, 0))));

	k = 0;
	for (int i = 0; i < dimset; i++) {
		for (int j = 0; j < dimset; j++) {
			complex<double> fourEnt(out[k][0] / dimset, out[k][1] / dimset); //Normalize results
			FourData.at(i).at(j).at(0) = fourEnt;
            if (pow(generateK(i, dimset, beam1.getK(), xKappa),2) + pow(generateK(j, dimset, beam1.getK(), yKappa),2) >= 1.0) //Catch complex wave vectors in z
				FourData.at(i).at(j).at(0) = (0,0); //Removes evanescence
			k++;
		}
	}

	//Recenter the Fourier data via 2D rotation
	
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

	//Allocate space for reflected beam E-field and calculate
	vector<vector<vector<complex<double> > > > eRTab(dimset, vector<vector<complex<double> > >(dimset, vector<complex<double> >(3, complex<double>(0, 0))));

	for (int j = 0; j < dimset; j++) {
		for (int i = 0; i < dimset; i++) {
			vector<double> kVec = vector<double> (0,0);
			//Calculate kappa vectors for each position in plane
			double kx = generateK(i + 1, dimset, beam1.getK(), xKappa); 
			double ky = generateK(j + 1, dimset, beam1.getK(), yKappa);

			kVec.push_back(kx);
			kVec.push_back(ky);

			//Catch complex kz's and set them to 0
			double kz = sqrt(1 - pow(kx, 2) - pow(ky, 2));
			if (isnan(kz)) kVec.push_back(0);
			else kVec.push_back(kz);
			eRTab.at(i).at(j) = eRBase(fVec, THETA, kVec);
		}
	}

	//cout << "Checkpoint: ETilde generated. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl; //Uncomment for timing studies

	vector<vector<vector<complex<double> > > > ERTab(dimset, vector<vector<complex<double> > >(dimset, vector<complex<double> >(3, complex<double>(0, 0))));
    int nanc = 0;
	
	//Mulitply by Fourier data from input beam to determine reflected beam in Fourier space
    for (int i = 0; i < dimset; i++) {
    	for (int j = 0; j < dimset; j++) {
			complex<double> fourPoint(FourData.at(i).at(j).at(0).real(), -FourData.at(i).at(j).at(0).imag());
			//Imaginary component gets flipped in FFTW
			
    		ERTab.at(i).at(j).at(0) = eRTab.at(i).at(j).at(0) * fourPoint;
    		ERTab.at(i).at(j).at(1) = eRTab.at(i).at(j).at(1) * fourPoint;
    		ERTab.at(i).at(j).at(2) = eRTab.at(i).at(j).at(2) * fourPoint;
			
			//Catch NaNs and set the corresponding points to 0
			if (isnan(ERTab.at(i).at(j).at(0).real()) || isnan(ERTab.at(i).at(j).at(1).real()) || isnan(ERTab.at(i).at(j).at(2).real()) || isnan(ERTab.at(i).at(j).at(0).imag()) || isnan(ERTab.at(i).at(j).at(1).imag()) || isnan(ERTab.at(i).at(j).at(2).imag())){
            cout << ++nanc << " NAN found: " << i << ", " << j << endl;
    		ERTab.at(i).at(j).at(0) = 0;
    		ERTab.at(i).at(j).at(1) = 0;
    		ERTab.at(i).at(j).at(2) = 0;
            }
                
			}
    }

	//cout << "Checkpoint: Ready for IFFT. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl; //Uncomment for timing studies

	in3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);
	out3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);
	
	//Fill IFT arrays
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

	//Use 3D transform to preserve extra dimensions we picked up during beam reconstruction. First two
	//dimensions are the same, but the third dimension is the number of components in each point (3)
	fftw_plan h3 = fftw_plan_dft_3d(ERTab.size(), ERTab.at(0).size(), 3, in3, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(h3);

	vector<vector<vector<complex<double> > > > outBeam(ERTab.size(), vector<vector<complex<double> > >(ERTab.at(0).size(), vector<complex<double> > (3, complex<double> (0,0))));
    vector<vector<vector<double> > > REoutBeam(ERTab.size(), vector<vector<double> > (ERTab.at(0).size(), vector<double> (3,0)));
	vector<vector<vector<double> > > outBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));
	vector<vector<vector<double> > > OGBeamMag(ERTab.size(), vector<vector<double> >(ERTab.at(0).size(), vector<double>(1, 0)));

	//cout << "After IFT:" << endl; //Uncomment for timing studies

	k = 0;
	for (int i = 0; i < ERTab.size(); i++) {
		for (int j = 0; j < ERTab.at(0).size(); j++) {
				vector<complex<double> > fourEnt;
				//Read the IFT matrix into a vector, normalize results, and store it
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
		}
	}

	//cout << "Checkpoint: Heavy lifting complete. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl; //Uncomment for timing studies

	
	//Cleanup 
	fftw_free(in);
	fftw_free(out);
	fftw_free(in3);
	fftw_free(out3);

	fftw_destroy_plan(h3);
    fftw_destroy_plan(g);

    //Calculate magnitude of both beams for comparision
	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		outBeamMag.at(i).at(j).at(0) = sqrt(real(outBeam.at(i).at(j).at(0) * conj(outBeam.at(i).at(j).at(0))) + real(outBeam.at(i).at(j).at(1) * conj(outBeam.at(i).at(j).at(1))) + real(outBeam.at(i).at(j).at(2) * conj(outBeam.at(i).at(j).at(2))));

	for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++) 
		OGBeamMag.at(i).at(j).at(0) = sqrt(pow(beam1.realEAt(i,j,0),2) + pow(beam1.imagEAt(i, j, 0),2));

    //Calculate centroid Shifts
	double nXrp1 = 0, nYrp1 = 0, denom = 0;
	
	for (int i = 0; i < outBeamMag.size(); i++) {
		for (int j = 0; j < outBeamMag.at(0).size(); j++) {
				nXrp1 += (j + 1) * pow(outBeamMag.at(j).at(i).at(0), 2);
				nYrp1 += (i + 1) * pow(outBeamMag.at(j).at(i).at(0), 2);
				denom += pow(outBeamMag.at(i).at(j).at(0), 2);
		}
	}

	double nXr = nXrp1 / denom;
	double nYr = nYrp1 / denom;
    double xShift = 2*400000/200*(nXr-(dimset+1)/2)/(2*PI/(632.8*pow(10, -9)))*pow(10,6);

    //Compare calculated shifts to analytical result
    complex<double> ARshift1 = (4*pow(NVAL,2)*sin(THETA*PI/180))/(beam1.getK()*(-1+pow(NVAL,2)+(1+pow(NVAL,2))*cos(2*THETA*PI/180))*sqrt(complex<double>(-1*pow(NVAL,2)+pow(sin(THETA*PI/180),2),0)));
    complex<double> ARshift2 = -2*sqrt(2)*sin(THETA*PI/180)/(beam1.getK()*sqrt(complex<double>(1-2*pow(NVAL,2)-cos(2*THETA*PI/180),0)));

	//cout << "Total time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;
	//Uncomment for timing studies

    return xShift;
}


int main(int argc, char** argv){
    int reso = 35;
    ofstream fout;
    vector<vector<double> > xShifts(reso+1, vector<double>(2,0));
    fout.open("GH_shifts.tsv");
    if(fout.fail()){ 
        cerr << "fout failed";
        exit(-1);
    }

    for (int i = 0; i <= reso; i++) {
		//Perform single-interface analysis across multiple angles and write to file
        double theta = 10+(70/(reso)*i);
        xShifts.at(i).at(0) = theta;
        xShifts.at(i).at(1) = shift(theta);
        fout << theta << "\t" << xShifts.at(i).at(1) << endl;
        cout << i << endl;
	}

    fout.close();
    return 0;
}
