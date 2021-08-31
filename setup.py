import numpy

'''
    setup.py file for spammodule.c

    Calling
    $python setup.py build_ext --inplace
    will build the extension library in the current file.

    Calling
    $python setup.py build
    will build a file that looks like ./build/lib*, where
    lib* is a file that begins with lib. The library will
    be in this file and end with a C library extension,
    such as .so

    Calling
    $python setup.py install
    will install the module in your site-packages file.

    See the distutils section of
    'Extending and Embedding the Python Interpreter'
    at docs.python.org for more information.
'''


from distutils.core import setup, Extension

canadafire_mod = Extension('canadafire', sources=['canadafire.c'],
    include_dirs=[numpy.get_include()])
#                        include_dirs=['/usr/local/lib'])

setup(name = 'canadafire',
        version='1.0',
        description='Canada fire weather XYZ',
        ext_modules = [canadafire_mod])
