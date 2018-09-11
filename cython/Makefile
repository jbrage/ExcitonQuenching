python=/usr/bin/python3

all:
	${python} setup.py build_ext --inplace 
	cython -a evolveDensitiesCython.pyx

clean:
	rm *.c *.so *.html
	rm -r build


