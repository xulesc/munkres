from distutils.core import setup, Extension

module1 = Extension('munkres', sources = ['munkres.c'], 
  extra_compile_args=['-std=c99', '-O3', '-shared', '-fPIC', '-ffast-math', '-lm', '-D __PYMOD__']) #-I/usr/include/python2.7/ -lpython2.7 -lm'])

setup (name = 'munkres',
  version = '0.1',
  description = 'kuhn-munkres algorithm for the Assignment Problem',
  license = "GPLv2",
  author = "Anuj Sharma",
  author_email = "anuj dot sharma80 at gmail dot com",
  ext_modules = [module1],
  classifiers = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GPLv2',
    'Operating System :: OS Independent',
    'Programming Language :: Python + C',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Software Development :: Libraries :: Python Modules'
  ]
)