# ExcitonQuenching
Calculation of ionization quenching correction factors in plastic scintillators.

Further details can be found in _Christensen JB and Andersen CE, Phys. Med. Biol. (2018)_  
https://doi.org/10.1088/1361-6560/aadf2d

_OBS: To be updated_  
Please contact jeppebrage@gmail.com if intended to use

## Installation
Linux (Tested on Ubuntu 16.04)

This version was tested with

* python2 and python3
* Cython 0.27 or newer 
* GNU Make v. 4.1

Install e.g. with

```
sudo apt-get install make
sudo apt-get install cython
```
## Run an example

Compile the cython code used to solve the Blanc equation
```
cd ExcitonQuenching/cython && make
```
which creates the shared object evolveDensitiesCython.so.
```
cd ../ && python main.py
```
will run a test version




