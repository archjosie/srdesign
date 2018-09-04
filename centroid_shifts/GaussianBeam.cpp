#include "GaussianBeam.h"

using namespace std;

GaussianBeam::GaussianBeam(double w0, double k, unsigned int p, unsigned int l) {
    this->w0 = w0;
    this->k = k;
    this->p = p;
    this->l = l;

    this->lambda = 2 * PI / k;
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
    // General vals that can be turned into real values
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

    // Populate yVals
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


    for (int i = 0; i < xVals.size(); ++i) {
        for (int j = 0; j < yVals.size(); ++j) {
            coordinate r = distance(xVals.at(i), yVals.at(j));

            // Calculate the real and imag arguments of exponential
            double reArg = -pow(k, 2)*pow(r, 2) / (pow(w0*k, 2));
            double imArg = 0;

            // get sign xval
            int sign = (xVals.at(i) > 0) - (xVals.at(i) < 0);
            if (abs(yVals.at(j)) < PI*1e-10) imArg = -l * (PI - sign * PI / 2);
            else imArg = -(l*atan(xVals.at(i) / yVals.at(j)));

            // Calculate beam intensity at all points in lattice
            intensity phasorOut;
            phasorOut = pow(2, l / 2)/ (w0 * k) * exp(reArg) * k * pow(k / (w0 * k * sqrt(1 / r)), l) * laguerre(p, l, 2 * pow(r, 2) * pow(k, 2) / (pow(w0*k, 2)));

            intensityReal realField = real(phasorOut);
            intensityImag imagField = imag(phasorOut);

            ReEField.at(i).at(j).at(0) = realField;
            ImEField.at(i).at(j).at(0) = imagField;
        }
    }
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
