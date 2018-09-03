#include "GaussianBeam.h"

using namespace std;

const static double NVAL = 0.659283;
const static double THETA = 60.0;

// Declarations
// ============

// Fourier Tranform Beam
// ---------------------
//beam<intensity> fourierTransformBeam(GaussianBeam beam1){
//    return
//}

// Helper functions
// ----------------

intensity rTE(double n, double theta, vector<double> kvec);
intensity rTM(double n, double theta, vector<double> kvec);
double generateK(int index, int dimsize, int k, int xmax);
double findMax(vector<double> vals);
vector<intensity> eRBase (vector<intensity> f, double theta, vector<double>  REkVecs);
double chop(double num);
int sign(double val);

// Main
// ====

int main(int argc, char** argv){
    time_t start = clock();
    cout << "Starting Calculations" << endl;

    // Initialize beam and relevent parameters
    double k0 = 1;
    GaussianBeam beam1(20000 / k0,2 * PI, 0, 0);
    beam1.calculateGaussData();
    int dimset = beam1.getDims();

    // Assuming horizontal polarization. According to Centroid Shifts paper, f={1,0,0}
    vector<intensity> fVec(3, intensity(0, 0));
    fVec.at(0) = 1;
    fVec.at(1) = complex<double>(0,-1);
    fVec.at(2) = 0;

    // Define the max xKappa and yKappa values
    double xKappa = findMax(beam1.getXVals());
    double yKappa = findMax(beam1.getYVals());

    // Generate the input and output vectors for fftw
    fftw_complex *in, *out, *in3, *out3;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimset * dimset);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimset * dimset);

    // Create plane for forward transform
    fftw_plan g = fftw_plan_dft_2d(dimset, dimset, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    int k = 0;
    for (int i = 0; i < dimset; i++) {
        for (int j = 0; j < dimset; j++) {
            in[k][0] = beam1.realEAt(i, j, 0);
            in[k][1] = beam1.imagEAt(i, j, 0);
            k++;
        }
    }

    cout << "Checkpoint: Fourier array prepared. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;
    
    // Do forward Fourier transform. Results are stored in "out"
    fftw_execute(g);

    // Extract Fourier data into a complex valued beam intensity
    beam<intensity> FourData(dimset, vector2D<intensity>  (dimset, vector<intensity>(1, intensity(0, 0))));

    k = 0;
    for (int i = 0; i < dimset; i++) {
        for (int j = 0; j < dimset; j++) {
            intensity fourEnt(out[k][0] / dimset, out[k][1] / dimset);
            FourData.at(i).at(j).at(0) = fourEnt;
            // Remove evanescence
            if (pow(generateK(i, dimset, beam1.getK(), xKappa), 2) + pow(generateK(j, dimset, beam1.getK(), yKappa), 2) >= 1.0)
                FourData.at(i).at(j).at(0) = (0,0);             
            k++;
        }
    }

/*   // Rotate three dimensional
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
*/
    // Calculate scalar that represents collision with beam
    beam<intensity> eRTab(dimset, vector2D<intensity> (dimset, vector<intensity>(3, intensity(0, 0))));

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
        }
    }

    cout << "Checkpoint: ETilde generated. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

    beam<intensity> ERTab(dimset, vector2D<intensity> (dimset, vector<intensity>(3, intensity(0, 0))));
    int nanc = 0;

    for (int i = 0; i < dimset; i++) {
        for (int j = 0; j < dimset; j++) {
            intensity fourPoint(FourData.at(i).at(j).at(0).real(), -FourData.at(i).at(j).at(0).imag());

            ERTab.at(i).at(j).at(0) = eRTab.at(i).at(j).at(0) * fourPoint;
            ERTab.at(i).at(j).at(1) = eRTab.at(i).at(j).at(1) * fourPoint;
            ERTab.at(i).at(j).at(2) = eRTab.at(i).at(j).at(2) * fourPoint;
            }
    }

    cout << "Checkpoint: Ready for IFFT. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

    in3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);
    out3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ERTab.size() * ERTab.at(0).size() * 3);

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

    fftw_plan h3 = fftw_plan_dft_3d(ERTab.size(), ERTab.at(0).size(), 3, in3, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(h3);

    beam<intensity>  outBeam(ERTab.size(), vector2D<intensity> (ERTab.at(0).size(), vector<intensity> (3, intensity (0,0))));
    beam<coordinate> REoutBeam(ERTab.size(), vector2D<intensityReal> (ERTab.at(0).size(), vector<intensityReal> (3,0)));
    beam<coordinate> outBeamMag(ERTab.size(), vector2D<double>(ERTab.at(0).size(), vector<double>(1, 0)));
    beam<coordinate> OGBeamMag(ERTab.size(), vector2D<double>(ERTab.at(0).size(), vector<double>(1, 0)));

    k = 0;
    for (int i = 0; i < ERTab.size(); i++) {
        for (int j = 0; j < ERTab.at(0).size(); j++) {
                vector<intensity> fourEnt;
                intensity entry1(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
                fourEnt.push_back(entry1);
                k++;

                intensity entry2(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
                fourEnt.push_back(entry2);
                k++;

                intensity entry3(out3[k][0] / (ERTab.size()*sqrt(3)), out3[k][1] / (ERTab.size()*sqrt(3)));
                fourEnt.push_back(entry3);
                k++;

                outBeam.at(i).at(j) = fourEnt;
        }
    }

    cout << "Checkpoint: Heavy lifting complete. " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

    fftw_free(in);
    fftw_free(out);
    fftw_free(in3);
    fftw_free(out3);

    fftw_destroy_plan(h3);
    fftw_destroy_plan(g);

    // Calculate magnitude of both beams for comparision
    for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++)
        outBeamMag.at(i).at(j).at(0) = sqrt(real(outBeam.at(i).at(j).at(0) * conj(outBeam.at(i).at(j).at(0))) + real(outBeam.at(i).at(j).at(1) * conj(outBeam.at(i).at(j).at(1))) + real(outBeam.at(i).at(j).at(2) * conj(outBeam.at(i).at(j).at(2))));

    for (int i = 0; i < outBeam.size(); i++) for (int j = 0; j < outBeam.at(0).size(); j++)
        OGBeamMag.at(i).at(j).at(0) = sqrt(pow(beam1.realEAt(i,j,0),2) + pow(beam1.imagEAt(i, j, 0),2));

    // Calculate centroid Shifts
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

    // Compare calculated shifts to analytical result
    intensity ARshift1 = (4*pow(NVAL,2)*sin(THETA*PI/180))/(beam1.getK()*(-1+pow(NVAL,2)+(1+pow(NVAL,2))*cos(2*THETA*PI/180))*sqrt(intensity(-1*pow(NVAL,2)+pow(sin(THETA*PI/180),2),0)));
    intensity ARshift2 = -2*sqrt(2)*sin(THETA*PI/180)/(beam1.getK()*sqrt(intensity(1-2*pow(NVAL,2)-cos(2*THETA*PI/180),0)));

    cout << "Calculated" << endl;
    cout << "(" << nXr << "," << nYr << ")" << endl;

    cout << "Total time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." << endl;

    return 0;
}

// Helper functions
// ----------------

intensity rTE(double n, double theta, vector<double> kvec) {
    theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
    double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);

    double beta = acos(theDot);
    intensity arg(pow(n, 2) - pow(sin(beta), 2), 0);
    intensity ans = (cos(beta) - sqrt(arg)) / (cos(beta) + sqrt(arg));

    return ans;
}

intensity rTM(double n, double theta, vector<double> kvec) {
    theta *= PI / 180;
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));
    double theDot = nvec.at(0)*kvec.at(0) + nvec.at(1)*kvec.at(1) + nvec.at(2)*kvec.at(2);

    double beta = acos(theDot);
    intensity arg(pow(n, 2) - pow(sin(beta), 2), 0);
    intensity ans = (pow(n,2)*cos(beta)-sqrt(arg))/(pow(n,2)*cos(beta)+sqrt(arg));


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

vector<intensity> eRBase (vector<intensity> f, double theta, vector<double> REkVecs) {
    vector<intensity> theVec(3, intensity (0,0));
    vector<intensity> kVecs(0, intensity (0,0));
    kVecs.push_back(intensity(REkVecs.at(0),0));
    kVecs.push_back(intensity(REkVecs.at(1),0));
    kVecs.push_back(intensity(REkVecs.at(2),0));
    intensity refTM = rTM(NVAL, theta, REkVecs);
    intensity refTE = rTE(NVAL, theta, REkVecs);

    theVec.at(0) = f.at(0)*refTM - f.at(1)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
    theVec.at(1) = f.at(1)*refTE + f.at(0)*kVecs.at(1)*(1 / tan(theta*PI/180))*(refTM + refTE);
    theVec.at(2) = -f.at(0)*refTM*kVecs.at(0) - f.at(1)*refTE*kVecs.at(1);

    return theVec;
}

double chop(double num) {
    if (abs(num) < 1e-12) return 0;
    return num;
}
