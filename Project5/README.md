# FYS4150 – Project 5

Folder containing code files used to solve Project 5 as well as the resulting
figures and text files.
----------------------------------------

Simulate2DSE.cpp/Simulate2DSE.hpp
--------
Class for simulate the 2-Dimensional Schrödinger Equation (2DSE) in a box, given
different parameters.

main.cpp
--------
C++ file that solves Project5 exercises using the Simulate2DSE class, writes files to directory "text_files"

plot.py
-------
Python file that plots results of exercises, saves pdf files to directory "figures"


Building and running the code:
-------------------------------
- main.cpp
     - Build: g++ *.cpp -o main.exe -O2 -std=c++11 -larmadillo
     - Run: ./main.exe
     - Build+run: make
- plot.py
     - Run: python plot.py
