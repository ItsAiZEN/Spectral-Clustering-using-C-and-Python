from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['spkmeansmodule.c'])
setup(name='mykmeanssp',
      version='1.0',
      description='A module that performs spectral clustering and returns wam, ddg, gl, jacobi matrices',
      ext_modules=[module])
