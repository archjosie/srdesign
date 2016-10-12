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
    // Set n value
    double n = 1.5;

    //Open root graphics
    TApplication theApp("App", &argc, argv);
    //Define canvas w/ dimensions
    TCanvas * c0 = new TCanvas("c1","multigraph L3",200,10,700,500);
    //Colors
    c0->SetFrameFillColor(30);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *gr1 = new TGraph(); gr1->SetLineColor(kBlue);
    TGraph *gr2 = new TGraph(); gr2->SetLineColor(kRed);
    TGraph *gr3 = new TGraph(); gr3->SetLineColor(kGreen);
    TGraph *gr4 = new TGraph(); gr4->SetLineColor(kOrange);
    //Make points for varied n values
    Double_t dx = 1;
    Double_t x  = 0;
    for (int i=0; i<=90; i++) {
       x = x+dx;
       gr1->SetPoint(i,x,transcofTE(n,x));
       gr2->SetPoint(i,x,transcofTE(n+.5,x));
       gr3->SetPoint(i,x,transcofTE(n+1.,x));
       gr4->SetPoint(i,x,transcofTE(n+1.5,x));
    }
    //Label Graph
    mg->Add(gr4); gr4->SetTitle("n=3.0"); gr4->SetLineWidth(3);gr4->SetFillStyle(0);
    mg->Add(gr3); gr3->SetTitle("n=2.5"); gr3->SetLineWidth(3);gr3->SetFillStyle(0);
    mg->Add(gr2); gr2->SetTitle("n=2.0"); gr2->SetLineWidth(3);gr2->SetFillStyle(0);
    mg->Add(gr1); gr1->SetTitle("n=1.5"); gr1->SetLineWidth(3);gr1->SetFillStyle(0);
    mg->SetTitle("Transverse Electric Coefficient of Transmission");
    mg->Draw("a");
    mg->GetXaxis()->SetTitle("Incident Angle (Degrees)"); 
    mg->GetYaxis()->SetTitle("Coefficient Value"); 
    mg->GetXaxis()->CenterTitle(); 
    mg->GetYaxis()->CenterTitle(); 
    //Create Legend
    TLegend *leg = new TLegend(.85, 0.7, .9, 0.9);
    leg->SetFillColor(0);
    leg->SetTextSize(0.02);
    leg->AddEntry(gr1,"n=1.5", "l");
    leg->AddEntry(gr2,"n=2.0", "l");
    leg->AddEntry(gr3,"n=2.5", "l");
    leg->AddEntry(gr4,"n=3.0", "l");
    leg->Draw();
    //Output PDF
    c0->Print("diffnplot.pdf","pdf");
    theApp.Run();
    return 0;
}
