# ExcitonQuenching
Calculation of ionization quenching correction factors in organic plastic scintillators.

Further details can be found in  
Christensen JB and Andersen CE (2018) _Relating ionization quenching in organic plastic scintillators to basic material properties by modelling excitation density transport and amorphous track structure during proton irradiation_, Phys. Med. Biol.      
https://doi.org/10.1088/1361-6560/aadf2d

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
cd ExcitonQuenching/EQ_cythonized_PDE && make
```
which creates the shared object evolveDensitiesCython.so.
```
cd ../ && python3 example.py
```
will execute an example




