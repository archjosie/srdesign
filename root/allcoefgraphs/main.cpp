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
#include <TLegend.h>
#include <TPad.h>
#include <TAxis.h>

#define PI 3.14159265

using namespace std;

// Set up equations for each coefficient
double refcofTE(double n,double thetaI){
    double arg, thetaTrad, thetaIrad, refcofTE; 
    thetaIrad = thetaI*PI/180;
    arg = 1/n * sin(thetaIrad);
    thetaTrad = asin(arg);
    refcofTE = abs((cos(thetaIrad) - n*cos(thetaTrad))/(cos(thetaIrad) + n*cos(thetaTrad))); 
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
    refcofTM = abs((cos(thetaTrad) - n*cos(thetaIrad))/(cos(thetaTrad) + n*cos(thetaIrad)));
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
    double n = 1.5;

    TApplication theApp("App", &argc, argv);
    TCanvas * c0 = new TCanvas("c1","multigraph L3",200,10,700,500);
    c0->SetFrameFillColor(30);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *gr1 = new TGraph(); gr1->SetLineColor(kBlue);
    TGraph *gr2 = new TGraph(); gr2->SetLineColor(kRed);
    TGraph *gr3 = new TGraph(); gr3->SetLineColor(kGreen);
    TGraph *gr4 = new TGraph(); gr4->SetLineColor(kOrange);
    Double_t dx = 1;
    Double_t x  = 0;
    for (int i=0; i<=90; i++) {
       x = x+dx;
       gr1->SetPoint(i,x,refcofTE(n,x));
       gr2->SetPoint(i,x,transcofTE(n,x));
       gr3->SetPoint(i,x,refcofTM(n,x));
       gr4->SetPoint(i,x,transcofTM(n,x));
    }
    mg->Add(gr4); gr4->SetLineWidth(3);gr4->SetFillStyle(0);
    mg->Add(gr3); gr3->SetLineWidth(3);gr3->SetFillStyle(0);
    mg->Add(gr2); gr2->SetLineWidth(3);gr2->SetFillStyle(0);
    mg->Add(gr1); gr1->SetLineWidth(3);gr1->SetFillStyle(0);
    mg->SetTitle("TE and TM Reflection and Transmission Coefficients");
    mg->Draw("a");
    mg->GetXaxis()->SetTitle("Incident Angle (Degrees)"); 
    mg->GetYaxis()->SetTitle("Coefficient Value"); 
    mg->GetXaxis()->CenterTitle(); 
    mg->GetYaxis()->CenterTitle(); 
    TLegend *leg = new TLegend(.55, 0.7, .75, 0.9);
    leg->SetFillColor(0);
    leg->SetTextSize(0.02);
    leg->AddEntry(gr1,"TE Reflection Coefficient", "l");
    leg->AddEntry(gr2,"TE Transmission Coeffecient", "l");
    leg->AddEntry(gr3,"TM Reflection Coefficient", "l");
    leg->AddEntry(gr4,"TM Transmission Coeffecient" , "l");
    leg->Draw();
    c0->Print("allcoeffsplot.pdf","pdf");
    theApp.Run();

    return 0;
}
