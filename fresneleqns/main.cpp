#include <iostream>
#include <math.h>
#include <cmath>
#include <assert.h> //Enables use of assert function (VERY useful for debugs)

#define PI 3.14159265

using namespace std;

//struct coeff {
//    double ref;
//    double trans;
//};

double snell(double thetaI, double n){
    double arg, thetaT, radtheta, thetaTdeg; 
	//Same as what you had, but on one line. Doubles default to 0.0 when no value is given
	
    radtheta = thetaI*PI/180;
    arg = 1/n * sin(radtheta);
    thetaT = asin(arg);
    thetaTdeg = asin(arg) *180/PI;
    cout << "Angle of Transmission " << thetaTdeg << endl;
    return thetaT;
}

void coeffs(double thetaT, double thetaI, double n){ 
//We aren't returning anything on this function, so the return statement should be void (for now)

    double refcofTE, refcofTM, transcofTE, transcofTM, refinTE, refinTM, transinTE, transinTM;
	
	//thetaI gets passed in degrees, so we need to fix that
	double thetaIrad = thetaI*PI/180;
	
	//The abs functions really don't do much, but I find positive numbers to be more intuitive to users
    refcofTE = abs((cos(thetaIrad) - n*cos(thetaT))/(cos(thetaIrad) + n*cos(thetaT))); 
    refinTE = pow (refcofTE,2);
    transcofTE = abs(2*cos(thetaIrad) /(cos(thetaIrad) + n*cos(thetaT)));
    transinTE= pow (transcofTE,2) * (n*cos(thetaT)/cos(thetaIrad)); //Snell factor included in calculation
    refcofTM = abs((cos(thetaT) - n*cos(thetaIrad))/(cos(thetaT) + n*cos(thetaIrad)));
    refinTM = pow (refcofTM,2);
    transcofTM = abs(2*cos(thetaIrad) /(cos(thetaT) + n*cos(thetaIrad)));
    transinTM = pow (transcofTM,2) * (n*cos(thetaT)/cos(thetaIrad)); //Snell factor included in calculation
	
	cout << "Transverse Electric Reflection Coefficient " << refcofTE << endl;
    cout << "Transverse Electric Transmission Coefficient " << transcofTE << endl;	
    cout << "Transverse Magnetic Reflection Coefficient " << refcofTM << endl;
    cout << "Transverse Magnetic Transmission Coefficient " << transcofTM << endl;
    cout << "Transverse Electric Reflection Intensity " << refinTE << endl;
    cout << "Transverse Electric Transmission Intensity " << transinTE << endl;	
    cout << "Transverse Magnetic Reflection Intensity " << refinTM << endl;
    cout << "Transverse Magnetic Transmission Intensity " << transinTM << endl;
}

int main(){
    double n, thetaI;
	//Other two variables were local variables we never use, so no need to declare them
	
    cout << "Value of n" << endl;
    cin >> n;
    cout << "Incident Angle" << endl;
    cin >> thetaI;
	
    double thetaT = snell(thetaI,n); //Declare thetaT and immediately assign it a value
    coeffs(thetaT, thetaI, n);
    return 0;
}
