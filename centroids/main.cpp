#include "GaussianBeam.h"

const static double PI = 3.14159265;
using namespace std;

int main(int argc, char** argv){
    GaussianBeam beam1(1.0,.01,1,0);
    beam1.calculateGaussData();
    //If we're only worrying about the interface we need to pick option 2
    //Flatten to 2D at a specific z    
    beam1.rootGraph(argc, argv, 0);
    return 0;
}
