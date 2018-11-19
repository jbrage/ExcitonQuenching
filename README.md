# ExcitonQuenching
Calculation of ionization quenching correction factors in organic plastic scintillators.

Details and background material given in  
Christensen JB and Andersen CE (2018) _Phys. Med. Biol._ __63__ 195010  
DOI: https://doi.org/10.1088/1361-6560/aadf2d

## Installation
Linux (Tested on Ubuntu 16.04)

This version was tested with

* python3
* Cython 0.27 or newer 
* GNU Make v. 4.1

Install e.g. with

```
sudo apt-get install make
sudo apt-get install cython
```

(Download e.g. with ```git clone https://github.com/jbrage/ExcitonQuenching```)

Compile the cython code used to solve the Blanc equation for a given ion and scintillator
```
cd ExcitonQuenching/EQ_cythonized_PDE && make
```
which creates the shared object evolveDensitiesCython.so.

## Run an example
Run an example with the BCF-12 scintillator exposed to a helium ion at several energies:  
```
python3 example.py
```




