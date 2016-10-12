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

#define PI 3.14159265

using namespace std;

void graph(char* function){
    // Define Plot and store in Draw?
    TF1 *fa1 = new TF1("fa1",function,0,90);
    fa1->Draw();
}

char* myf(double n){
    //construct the function string that will be input to the graph function
    char buffer [50], sbuffer [1000];
    sprintf(buffer,"asin(1/ %4.2f * sin(x*pi/180))", n);
    sprintf(sbuffer,"(cos(x*pi/180) - %4.2f*cos(%s))/(cos(x*pi/180) + %4.2f*cos(%s))", n, buffer, n, buffer); 

//    std::ostringstream oss;
//    oss << sbuffer;
//    string funcstr = oss.str();

    return sbuffer;
}

int main(int argc, char *argv[]){
    double n, thetaI;
    cout << "Value of n" << endl;
    cin >> n;

    char buffer [50], sbuffer [1000];
    sprintf(buffer,"asin(1/ %4.2f * sin(x*pi/180))", n);
    sprintf(sbuffer,"(cos(x*pi/180) - %4.2f*cos(%s))/(cos(x*pi/180) + %4.2f*cos(%s))", n, buffer, n, buffer); 

    TApplication theApp("App", &argc, argv);
    graph(sbuffer);
    theApp.Run();

    return 0;
}
