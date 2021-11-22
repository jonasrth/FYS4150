# FYS4150 - Project 4
Repository for Project 4

Contains files needed to solve Project 4

main.cpp
--------
C++ file that solves Project4 exercises, writes to terminal and writes files to directory "text_files"

plot.py
-------
Python file that plots results of exercises, saves pdf files to directory "figures"


Building and running the code:
-------------------------------
- main.cpp
     Build: g++ main.cpp RandomFlipper.cpp -Xpreprocessor -fopenmp -o main.exe -O2 -std=c++11 -larmadillo -lomp
     Run: ./main.exe
     Build+run: make
- plot.py
     Run: python plot.py
