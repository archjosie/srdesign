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

    double calculateWaist(double z);
    double calculateRadCurv(double z);
    double calculateRayleigh();
    double calculateGouy(double z);
    double calculateHermite(double x, unsigned int m);

public:
    // paramaterized constructor
    GaussianBeam(double w0, double lambda, unsigned int m, unsigned int n);
    // default constructor
    GaussianBeam() : GaussianBeam(1, 20, 0, 0) {};

    void calculateGaussData();
    void fourierTran(beam<intensityReal> realPart, beam<intensityImag> imagPart, vector<double> xVals, vector<double> yVals, ofstream &fout);
    static double distance(double x, double y);
    void setRealE(beam<intensityReal> realE);
    void setImE(beam<intensityImag> imE);
    void setRealFour(beam<intensityReal> realFour);
    void setImFour(beam<intensityImag> imFour);
    beam<intensityReal> getRealE();
    beam<intensityImag> getImE();
    vector<coordinate> getXVals();
    vector<coordinate> getYVals();
    double getK();
    double laguerre(unsigned int k, double alpha, double x);
    int getDims();
    double realEAt(int i, int j, int k);
    double imagEAt(int i, int j, int k);
};

#endif
