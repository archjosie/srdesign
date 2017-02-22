#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <assert.h> 

#define PI 3.14159265

using namespace std;

double rTE(double n, double theta, vector<double> kvec){
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (cos(beta)-sqrt(pow(n,2)-pow(sin(beta),2)))/(cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)));
    return coeff;
} 

double rTM(double n, double theta, vector<double> kvec){
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (pow(n,2)*cos(beta)-sqrt(pow(n,2)-pow(sin(beta),2)))/(pow(n,2)*cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)));
    return coeff;
} 

double tTE(double n, double theta, vector<double> kvec){
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (2*cos(beta))/(cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)));
    return coeff;
} 

double tTM(double n, double theta, vector<double> kvec){
    vector<double> nvec;
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (2*n*cos(beta))/(pow(n,2)*cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)));
    return coeff;
} 

int main(){
//    double n = 1.6592;
    double n = .6;
    double theta = 30*PI/180;
    vector<double> kvec;
    kvec.push_back(0.0);
    kvec.push_back(0.0);
    kvec.push_back(1.0);

    cout << rTE(n, theta, kvec) << endl;
    cout << rTM(n, theta, kvec) << endl;
    cout << tTE(n, theta, kvec) << endl;
    cout << tTM(n, theta, kvec) << endl;
	
    return 0;
}
