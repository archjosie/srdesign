#include "GaussianBeam.h"

const static double PI = 3.14159265;
using namespace std;

int main(int argc, char** argv){
    GaussianBeam beam1(1.0,.01,1,0);
    beam1.calculateGaussData();
    beam1.rootGraph(argc, argv, 0);
    return 0;
}
