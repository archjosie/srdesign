#include "GaussianBeam.h"

const static double PI = 3.14159265;



void Graph(int argc, char** argv, vector<vector<double> > theVec, vector<double> xVals, vector<double> yVals){
    //Open root graphics
    TApplication theApp("App", &argc, argv);
    gStyle->SetOptStat(0);
//    gStyle->SetPalette(82);
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    TH2F *hcontz = new TH2F("hcontz","Gaussian Beam Cross Section",40,-PI,PI,40,-PI,PI);
	for (int i = 0; i < xVals.size(); ++i) {
		for (int j = 0; j < yVals.size(); ++j) {
				hcontz->Fill(xVals.at(i),yVals.at(j),theVec.at(i).at(j));
		}
	}
//    hcontz->SetTitle("Cross Section of Gaussian Beam");
    hcontz->GetXaxis()->SetTitle("x"); 
    hcontz->GetYaxis()->SetTitle("y"); 
    hcontz->SetMarkerStyle(1);
    hcontz->GetXaxis()->CenterTitle(); 
    hcontz->GetYaxis()->CenterTitle(); 
    
    hcontz->Draw("CONTZ");
    // Output PDF
    c1->Print("GBplots.pdf","pdf");
    theApp.Run();
}

int main(int argc, char** argv){

    int sizeX = 2000 ;
    int sizeY = 2000 ;
    
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);

	fftw_plan g = fftw_plan_dft_2d(sizeX, sizeY, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	int k = 0;
	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < sizeY; j++) {
			in[k][0] = 0;//sin(5 * PI/sizeY*j);//+ sin(10*PI/sizeX*j);//+sin (.5* PI/sizeX*i)+sin (100*PI/sizeX*i)+sin (.25*PI/sizeX*i)+sin (.1*PI/sizeX*i)+sin (PI/sizeX*i); 
			in[k][1] = sin(PI/sizeY*j);//0;//sin (10*PI/sizeY*j)+sin (.5* PI/sizeX*j)+sin (100*PI/sizeX*j)+sin (.25*PI/sizeX*j)+sin (.1*PI/sizeX*j)+sin (PI/sizeX*j);
			k++;
		}
	}

	fftw_execute(g);

    vector<vector<double> >  outVec(sizeX, vector<double> (sizeY, 0));
    vector<double>  xVals(sizeX, 0);
    vector<double>  yVals(sizeY, 0);
	//vector<vector<vector<double> > > ImFour(sizeX, vector<vector<double> >(sizeY, vector<double>(1, 0)));

	k = 0;
	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < sizeY; j++) {
			outVec.at(i).at(j) = out[k][0];
            xVals.at(i) = sin (10*PI/sizeX*i);
            yVals.at(j) = sin (10*PI/sizeY*j);
    //		ImFour.at(i).at(j).at(0) = out[k][1];
			k++;
		}
	}

    Graph(argc, argv, outVec, xVals, yVals);

    return 0;
}
