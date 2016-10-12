#include <iostream>
#include <math.h>

#define PI 3.14159265

using namespace std;

int main(){
    double n1, n2, theta, rangle, arg, radtheta;
    string unit;

    cout << "Value of n1" << endl;
    cin >> n1;
    cout << "Value of n2" << endl;
    cin >> n2;
    cout << "Value of theta" << endl;
    cin >> theta;
    cout << "Units of theta (rad or deg)" << endl;
    cin >> unit;
   
    if(unit == "deg"){
        radtheta = theta*PI/180;
        arg = n1/n2* sin(radtheta);
        rangle = asin(arg) *180/PI;
        cout << "Angle of refraction " << rangle << endl;
    } else if(unit == "rad"){
        arg = n1/n2* sin(theta);
        rangle = asin(arg);
        cout << "Angle of refraction " << rangle << endl;
    } else{
        cout << "Units not recognized " << endl;
    }
	
    return 0;
}
