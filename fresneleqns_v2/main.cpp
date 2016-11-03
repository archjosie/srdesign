#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <assert.h> //Enables use of assert function (VERY useful for debugs)

#define PI 3.14159265

using namespace std;

double snell(double n, double thetaIrad){
    double arg, thetaTrad; 
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    return thetaTrad;
}

double refcofTE(double n,double thetaI){
    double thetaTrad, thetaIrad, refcofTE; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    refcofTE = abs((cos(thetaIrad) - n*cos(thetaTrad))/(cos(thetaIrad) + n*cos(thetaTrad))); 
    return refcofTE;
}

double refinTE(double n,double thetaI){
    double thetaTrad, thetaIrad, refcofTE, refinTE; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    refcofTE = abs((cos(thetaIrad) - n*cos(thetaTrad))/(cos(thetaIrad) + n*cos(thetaTrad))); 
    refinTE = pow(refcofTE,2);
    return refinTE;
}

double transcofTE(double n,double thetaI){
    double thetaTrad, thetaIrad, transcofTE; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    transcofTE = (2*cos(thetaIrad) /(cos(thetaIrad) + n*cos(thetaTrad)));
    return transcofTE;
}

double transinTE(double n,double thetaI){
    double thetaTrad, thetaIrad, transcofTE, transinTE; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    transcofTE = (2*cos(thetaIrad) /(cos(thetaIrad) + n*cos(thetaTrad)));
    transinTE= pow (transcofTE,2) * (n*cos(thetaTrad)/cos(thetaIrad)); //Snell factor included in calculation
    return transinTE;
}

double refcofTM(double n,double thetaI){
    double thetaTrad, thetaIrad, refcofTM; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    refcofTM = abs((cos(thetaTrad) - n*cos(thetaIrad))/(cos(thetaTrad) + n*cos(thetaIrad)));
    return refcofTM;
}

double refinTM(double n,double thetaI){
    double thetaTrad, thetaIrad, refcofTM, refinTM; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    refcofTM = abs((cos(thetaTrad) - n*cos(thetaIrad))/(cos(thetaTrad) + n*cos(thetaIrad)));
    refinTM = pow (refcofTM,2);
    return refinTM;
}


double transcofTM(double n,double thetaI){
    double thetaTrad, thetaIrad, transcofTM; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    transcofTM = (2*cos(thetaIrad) /(cos(thetaTrad) + n*cos(thetaIrad)));
    return transcofTM;
}

double transinTM(double n,double thetaI){
    double thetaTrad, thetaIrad, transcofTM, transinTM; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    transcofTM = (2*cos(thetaIrad) /(cos(thetaTrad) + n*cos(thetaIrad)));
    transinTM = pow (transcofTM,2) * (n*cos(thetaTrad)/cos(thetaIrad)); //Snell factor included in calculation
    return transinTM;
}

void wavenums(double k, double n, double thetaI, vector<double> &kin, vector<double> &kref, vector<double> &ktran){
    double thetaIrad, thetaTrad; 
    thetaIrad = thetaI*PI/180;
    thetaTrad = snell(n, thetaIrad);
    
    kin.push_back(sin(thetaIrad));
    kin.push_back(0);
    kin.push_back(cos(thetaIrad));
    
    kref.push_back(-sin(thetaIrad));
    kref.push_back(0);
    kref.push_back(cos(thetaIrad));
    
    ktran.push_back((1/n)*sin(thetaTrad));
    ktran.push_back(0);
    ktran.push_back((1/n)*cos(thetaTrad));
}

void coeffs(double n, double thetaI){ 
    vector<double> kin, kref, ktran;
    wavenums(n, thetaI, kin, kref, ktran);
	cout << "Transverse Electric Reflection Coefficient " <<    refcofTE(n,thetaI) << endl;
    cout << "Transverse Electric Transmission Coefficient " <<  transcofTE(n,thetaI) << endl;	
    cout << "Transverse Magnetic Reflection Coefficient " <<    refcofTM(n,thetaI) << endl;
    cout << "Transverse Magnetic Transmission Coefficient " <<  transcofTM(n,thetaI) << endl;
    cout << "Transverse Electric Reflection Intensity " <<      refinTE(n,thetaI) << endl;
    cout << "Transverse Electric Transmission Intensity " <<    transinTE (n,thetaI) << endl;	
    cout << "Transverse Magnetic Reflection Intensity " <<      refinTM(n,thetaI) << endl;
    cout << "Transverse Magnetic Transmission Intensity " <<    transinTM(n,thetaI) << endl;
    cout << "Incident Wave number: ("<< kin.at(0) << "," << kin.at(1) << "," << kin.at(2) << ")" << endl;
    cout << "Reflected Wave number: ("<< kref.at(0) << "," << kref.at(1) << "," << kref.at(2) << ")" << endl;
    cout << "Transmitted Wave number: ("<< ktran.at(0) << "," << ktran.at(1) << "," << ktran.at(2) << ")" << endl;
}

int main(){
    double n, thetaI, k;
	//Other two variables were local variables we never use, so no need to declare them
	
    cout << "Value of n" << endl;
    cin >> n;
    cout << "Incident Angle" << endl;
    cin >> thetaI;
    cout << "Wave Number" << endl;
    cin >> k;
	
    coeffs(n, thetaI);
    return 0;
}
