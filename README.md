# ExcitonQuenching
Calculation of ionization quenching correction factors in plastic scintillators


# Installation
Linux (Tested on Ubuntu 16.04)

This version was tested with

Python2.7
Cython 0.27 or newer 
GNU Make v. 4.1

sudo apt-get install make
sudo apt-get install python-pip
pip install Cython

# Compile and run an example

cd ExcitonQuenching/cython
make

(should create the share object evolveDensitiesCython.so)

cd ../
python2.7 main.py

(should run a test version)



