#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <assert.h> //Enables use of assert function (VERY useful for debugs)

#define PI 3.14159265

using namespace std;

double rTE(double n, double theta, vector<double> kvec){
    vector<double> nvec(3,0);
    nvec.push_back(-sin(theta));
    nvec.push_back(0);
    nvec.push_back(cos(theta));

    double beta= acos(nvec.at(0)*kvec.at(0)+nvec.at(1)*kvec.at(1)+nvec.at(2)*kvec.at(2));
    double coeff= (cos(beta)-sqrt(pow(n,2)-pow(sin(beta),2)))/(cos(beta)+sqrt(pow(n,2)-pow(sin(beta),2)))
} 


int main(){
    double n = 0.6592;
    double theta = 30*PI/180;
    vector<double> kvec(3,0);
    kvec.push_back(0);
    kvec.push_back(0);
    kvec.push_back(1);




	
    return 0;
}
