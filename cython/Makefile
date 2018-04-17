all:
	python2.7 setup.py build_ext --inplace 
	cython -a evolveDensitiesCython.pyx
#	python main.py 

clean:
	rm *.c *.so *.html
	rm -r build


