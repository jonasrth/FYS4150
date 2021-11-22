all: compile execute

compile:
	#g++ *.cpp -Xpreprocessor -fopenmp -o main.exe -O2 -std=c++11 -larmadillo -lomp
	g++ main.cpp RandomFlipper.cpp -Xpreprocessor -fopenmp -o main.exe -O2 -std=c++11 -larmadillo -lomp

execute:
	./main.exe
