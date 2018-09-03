#include "GaussianBeam.h"

using namespace std;

GaussianBeam::GaussianBeam(double w0, double lambda, unsigned int p, unsigned int l) {
    this->w0 = w0;
    this->lambda = lambda;
    this->p = p;
    this->l = l;

    this->k = 2 * PI / lambda;
    this->zR = calculateRayleigh();
}

double GaussianBeam::distance(double x, double y) {
    return sqrt(pow(x, 2) + pow(y, 2));
}

double GaussianBeam::calculateRadCurv(double z) {
    return z * (1 + pow(zR / z, 2));
}

double GaussianBeam::calculateGouy(double z) {
    return atan(z / zR) * (l + 2 * p);
}

double GaussianBeam::calculateWaist(double z) {
    return w0 * sqrt(1 + pow(z / zR, 2));
}

double GaussianBeam::calculateRayleigh() {
    return PI * pow(w0, 2) / lambda;
}

// Calculates the nth Hermite polynomial evaluated at x recursively
double GaussianBeam::calculateHermite(double x, unsigned int m) {
    // Base cases
    if (m == 0) return 1;
    if (m == 1) return 2 * x;

    return 2 * x * calculateHermite(x, m - 1) - 2 * (m - 1) * calculateHermite(x, m - 2);
}

void GaussianBeam::calculateGaussData() {
    dimset = 201;
    xMax = 400000;
    xMin = -xMax;
    xInt = 2 * xMax / (dimset - 1);
    yMax = xMax;
    yMin = -yMax;
    yInt = xInt;

    int xRange = (xMax - xMin) / xInt + 1;
    int yRange = (yMax - yMin) / yInt + 1;

    // Populate xVals
    for (int i = 0; i < xRange; i++) {
        double xCurr = xMin + i*xInt;
        xVals.push_back(xCurr);
    }

    // Populate rVals
    for (int i = 0; i < yRange; i++) {
        double yCurr = yMin + i*yInt;
        yVals.push_back(yCurr);
    }

    // Stores lower limit on distance from focus
    double theZ;
    theZ = .1;
    zVals.push_back(theZ);

    // Separate 3D vectors to hold real and imaginary parts of E-field
    beam<intensityReal> ReLocal(xVals.size(), vector2D<intensityReal>(yVals.size(), vector<intensityReal>(zVals.size(), 0)));

    ReEField = ReLocal;
    ImEField = ReLocal;

    double omega = 20 / k;

    for (int i = 0; i < xVals.size(); ++i) {
        for (int j = 0; j < yVals.size(); ++j) {
            coordinate r = distance(xVals.at(i), yVals.at(j));

            intensityReal reArg = -pow(k, 2)*pow(r, 2) / (pow(omega*k, 2));
            intensityImag imArg = 0;
            // get sign xval
            int sign = (xVals.at(i) > 0) - (xVals.at(i) < 0);
            if (abs(yVals.at(j)) < PI*1e-10) imArg = -l * (PI - sign * PI / 2);
            else imArg = -(l*atan(xVals.at(i) / yVals.at(j)));

            // Uncomment if we only want to consider real component of field
            intensity phasorOut;
            phasorOut = pow(2, l / 2)/ (omega * k) * exp(reArg) * k * pow(k / (omega * k * sqrt(1 / r)), l) * laguerre(p, l, 2 * pow(r, 2) * pow(k, 2) / (pow(omega*k, 2)));

            intensityReal realField = real(phasorOut);
            intensityImag imagField = imag(phasorOut);

            ReEField.at(i).at(j).at(0) = realField;
            ImEField.at(i).at(j).at(0) = imagField;
        }
    }
}

void GaussianBeam::fourierTran(beam<intensityReal> realPart, beam<intensityImag> imagPart, vector<coordinate> xVals, vector<coordinate> yVals, ofstream &fout) {
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

    vector2D<intensityReal> realFour(realPart.size(), vector<intensityReal>(realPart.at(0).size()));
    vector2D<intensityImag> imagFour(imagPart.size(), vector<intensityImag>(imagPart.at(0).size()));

    k = 0;
    for (int i = 0; i < realPart.size(); ++i) {
        for (int j = 0; j < realPart.at(0).size(); ++j) {
            realFour.at(i).at(j) = out[k][0];
            imagFour.at(i).at(j) = out[k][1];
            k++;
        }
    }

    // Pound sign forces gnuplot to ignore line
    // Write our 3D complex vector to file
    fout << "# x " << '\t' << "y " << '\t' << "Re(E(x,y)) " << '\t' << "Im(E(x,y)) " << '\t' << "Abs(E(x,y))" << endl;
    for (int i = 0; i < xVals.size(); ++i) {
        for (int j = 0; j < yVals.size(); ++j) {
                fout << xVals.at(i) << '\t' << yVals.at(j) << '\t' << realFour.at(i).at(j) << '\t' << imagFour.at(i).at(j)<< '\t' << distance(realFour.at(i).at(j), imagFour.at(i).at(j)) << endl;
            }
        }

    cout << "Data output complete!" << endl;
}

beam<intensityReal> GaussianBeam::getRealE() {
    return ReEField;
}

beam<intensityImag> GaussianBeam::getImE() {
    return ImEField;
}

vector<coordinate> GaussianBeam::getXVals() {
    return xVals;
}

vector<coordinate> GaussianBeam::getYVals() {
    return yVals;
}

double GaussianBeam::getK() {
    return k;
}

double GaussianBeam::laguerre(unsigned int k, double alpha, coordinate x) {
    if (k == 0) return 1;
    if (k == 1) return 1 + alpha - x;

    return ((2 * k + 1 + alpha - x) * laguerre(k - 1, alpha, x) - (k + alpha) * laguerre(k - 2, alpha, x)) / (k + 1);
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
