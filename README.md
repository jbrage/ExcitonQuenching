# ExcitonQuenching
Calculation of ionization quenching correction factors in plastic scintillators

_OBS: To be updated_  
Please contact jeppebrage@gmail.com if intended to use

## Installation
Linux (Tested on Ubuntu 16.04)

This version was tested with

* Python2.7
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
cd ../ && python2 main.py
```
will run a test version




