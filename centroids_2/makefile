default: main

main: 
	g++ -o main main.cpp GaussianBeam.cpp -lfftw3 -lm `root-config --cflags --libs`

clean:
	rm -f main
