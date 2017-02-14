#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char** argv){
    ifstream din;
    ofstream dout;
    din.open("ddpostprocess_1.out");
	
	if (din.fail()) {
		cerr << "DD post-process file not found." << endl;
		exit(1);
	}

	string dummy;
    for (int i = 0; i < 17; i++) getline(din, dummy); //Eliminate filler lines and labels

    vector<double> x,y,z,rx,ix,ry,iy,rz,iz; //Set up vectors to hold data
    while(!din.eof()){
		double tmp;
		din >> tmp;
        x.push_back(tmp); //Magnitude x

		din >> tmp;
		y.push_back(tmp); //Magnitude y

		din >> tmp;
		z.push_back(tmp); //Magnitude z
		
		din >> tmp;
		rx.push_back(tmp); //Real part of x

		din >> tmp;
		ix.push_back(tmp); //Real part of y

		din >> tmp;
		ry.push_back(tmp); //Real part of z

		din >> tmp;
		iy.push_back(tmp); //Imaginary part of x

		din >> tmp;
		rz.push_back(tmp); //Imaginary part of y

		din >> tmp;
		iz.push_back(tmp); //Imaginary part of z
    }

	for (int i = 0; i < 20; i++) cout << x.at(i) << endl; //Print 20 lines of x as test
   
	din.close();
  
    return 0;
}
