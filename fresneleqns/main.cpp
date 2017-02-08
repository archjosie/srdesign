#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <assert.h> 

#define PI 3.14159265

using namespace std;

double rTE(double n, double theta, vector<double> kvec){
    vector<double> nvec(3);
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (cos(beta)-sqrt(pow(n,2)-pow(sin(beta),2)))/(cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)));
    return coeff;
} 

double rTEb(double n, double theta, vector<double> kvec){
    vector<double> nvec(3);
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    return beta;
}

int main(){
//    double n = 1.6592;
    double n = .6;
    double theta = 30*PI/180;
    vector<double> kvec(3);
    kvec.push_back(0.0);
    kvec.push_back(0.0);
    kvec.push_back(1.0);

    cout << kvec.at(3) << endl;
    cout << rTE(n, theta, kvec) << endl;
    cout << rTEb(n, theta, kvec) << endl;

	
    return 0;
}
