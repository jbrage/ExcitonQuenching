from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = 'ExcitonQuenching',
        description = 'Quenching correction factors for organic scintillators exposed to ions',
        author = 'Jeppe Brage Christensen',
        author_email = 'jeppebrage@gmail.com',
        url = 'https://github.com/jbrage/ExcitonQuenching',
        ext_modules=cythonize("evolveDensitiesCython.pyx"),
)
