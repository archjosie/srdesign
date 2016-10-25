#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <assert.h> //Enables use of assert function (VERY useful for debugs)
#include <sstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TMath.h>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TAxis.h>

#define PI 3.14159265

using namespace std;

//Set up equations for each different coefficient
double refcofTE(double n,double thetaI){
    double arg, thetaTrad, thetaIrad, refcofTE; 
    thetaIrad = thetaI*PI/180;
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    refcofTE = ((cos(thetaIrad) - n*cos(thetaTrad))/(cos(thetaIrad) + n*cos(thetaTrad))); 
    return refcofTE;
}

double transcofTE(double n,double thetaI){
    double arg, thetaTrad, thetaIrad, transcofTE; 
    thetaIrad = thetaI*PI/180;
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    transcofTE = (2*cos(thetaIrad) /(cos(thetaIrad) + n*cos(thetaTrad)));
    return transcofTE;
}

double refcofTM(double n,double thetaI){
    double arg, thetaTrad, thetaIrad, refcofTM; 
    thetaIrad = thetaI*PI/180;
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    refcofTM = ((cos(thetaTrad) - n*cos(thetaIrad))/(cos(thetaTrad) + n*cos(thetaIrad)));
    return refcofTM;
}

double transcofTM(double n,double thetaI){
    double arg, thetaTrad, thetaIrad, transcofTM; 
    thetaIrad = thetaI*PI/180;
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    transcofTM = (2*cos(thetaIrad) /(cos(thetaTrad) + n*cos(thetaIrad)));
    return transcofTM;
}

int main(int argc, char *argv[]){
    //Open root graphics
    TApplication theApp("App", &argc, argv);
    //Define canvas w/ dimensions
    TCanvas *c = new TCanvas("c","Graph2D example",0,0,600,400);
    //Use T class doubles to make random points on the surface
    Double_t n, thetaI, z;
    int np = 3600;
    TGraph2D *dt = new TGraph2D();
    //Create range of n values
    TF1 *f2 = new TF1("f2","x",1,3);
    for (Int_t N=0; N<np; N++) {
        // theta is iternated in .1deg intervals and n is chosen randomly from our range
        thetaI = N/40;
        n = f2->GetRandom();
        z = transcofTM(n,thetaI);
        dt->SetPoint(N,thetaI,n,z);
    }
    gStyle->SetPalette(1);
    // Set title and Axis Labels
    dt->SetTitle("TM Transmission Coefficients at Varied Angles of Incident and Indicies of Refraction");
    dt->GetXaxis()->SetTitle("Incident Angle (Degrees)"); 
    dt->GetYaxis()->SetTitle("Coefficient Value"); 
    dt->SetMarkerStyle(21);
    dt->GetXaxis()->CenterTitle(); 
    dt->GetYaxis()->CenterTitle(); 
    dt->Draw("PCOL");
    // Output PDF
    c->Print("3dplots.pdf","pdf");
    theApp.Run();

    return 0;
}
