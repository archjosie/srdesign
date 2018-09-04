#ifndef GAUSSIANBEAM_H
#define GAUSSIANBEAM_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <fstream>
#include <vector>
#include <complex>

#include<fftw3.h>

using namespace std;

// type synonyms
template <typename T>
using beam = vector<vector<vector<T>>>;
template <typename T>
using vector2D = vector<vector<T>>;
using coordinate = double;
using intensity = complex<double>;
using intensityReal = double;
using intensityImag = double;

const static double PI = 3.14159265;

class GaussianBeam {

// Comment
// - I think some of the private functions here aren't defined in the .cpp,
//   should be deleted if so. Maybe some of the public too, I haven't checked

private:
    double w0, lambda, k, zR;
    unsigned int p, l;
    double dimset;
    double xMin, xMax, xInt;
    double yMin, yMax, yInt;
    double zMin, zMax, zInt;
    double tMin, tInt, tRange;

    beam<intensityReal> ReEField;
    beam<intensityImag> ImEField;
    vector<double> rVals, tVals, zVals, xVals, yVals;

    static double distance(double x, double y);
    double calculateRadCurv(double z);
    double calculateGouy(double z);
    double calculateWaist(double z);
    double calculateRayleigh();
    double calculateHermite(double x, unsigned int m);
    double laguerre(unsigned int k, double alpha, double x);

public:
    // parameterized constructor
    GaussianBeam(double w0, double lambda, unsigned int m, unsigned int n);
    // default constructor
    GaussianBeam() : GaussianBeam(1, 20, 0, 0) {};

    void calculateGaussData();
    vector<coordinate> getXVals();
    vector<coordinate> getYVals();
    double getK();
    int getDims();
    double getLen();
    double realEAt(int i, int j, int k);
    double imagEAt(int i, int j, int k);
};

#endif
