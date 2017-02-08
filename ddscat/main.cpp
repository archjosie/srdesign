#include "GaussianBeam.h"

int main(int argc, char** argv){
    ifstream din;
    ofstream dout;
    din.open("ddpostprocess_1.out");

    if (din.is_open()){
        for(int i=0; i<17; i++){
           getline(din, string line);
        }
        vector<double> x,y,z,rx,ix,ry,iy,rz,iz;
        while(!din.eof()){
            x.push_back()=stod(getline(din, string data, "  "));
            y.push_back()=stod(getline(din, string data, "  "));
            z.push_back()=stod(getline(din, string data, "  "));
            rx.push_back()=stod(getline(din, string data, "  "));
            ix.push_back()=stod(getline(din, string data, "  "));
            ry.push_back()=stod(getline(din, string data, "  "));
            iy.push_back()=stod(getline(din, string data, "  "));
            rz.push_back()=stod(getline(din, string data, "  "));
            iz.push_back()=stod(getline(din, string data, "  "));
        }
        for(int i=0; i<20; i++){
            cout << x.at(i) << endl;
        }
    }
  
    return 0;
}
