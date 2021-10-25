# FYS4150 - Project 3
Repository for Project 3

main.cpp
--------
C++ file that solves Project3 exercises, writes to terminal and writes files to directory "text_files"

plot.py
-------
Python file that plots results of exercises, saves pdf files to directory "figures"


Building and running the code:
-------------------------------
- main.cpp
     Build: g++ main.cpp *.cpp -o main.exe -std=c++11 -O2 -larmadillo
     Run: ./main.exe

     Build+run: make
- plot.py
     Run: python plot.py
