"""A setuptools based setup module.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension, Command
from os import system

# Numpy is always needed
try:
    import numpy
except ImportError:
    print("Install numpy first, required for compilation:")
    print("    pip3 install numpy")
    print("or")
    print("    pip3 install --user numpy")
    exit()

# Cython is useful
try:
    import Cython
except ImportError:
    print("No Cython found.")
    print(" Continuing, but installation of Cython is recommended")
    use_cython = False
else:
    use_cython = True

if use_cython:
    ext_modules = [Extension("TumOnc.critstat", ["TumOnc/critstat.pyx"])]
else:
    ext_modules = [Extension("TumOnc.critstat", ["TumOnc/critstat.c"])]


class MinusPublish(Command):
    description = 'publish to minus web server (internal use only)'
    user_options = []

    def initialize_options(self):
        return

    def finalize_options(self):
        return

    def run(self):
        system("rst2html.py README.rst >README.html")
        system("rsync -avP README.html dist/* minus:public_html/TumOnc/")


setup(
    name='TumOnc',
    version='0.0.1.dev5',
    url='http://www.inr.ac.ru/~fedor/TumOnc',
    author='Fedor Bezrukov',
    author_email='Fedor.Bezrukov@gmail.com',
    ext_modules=ext_modules,
    include_dirs=[numpy.get_include()],
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'scipy', 'matplotlib'],
    scripts=['TumOnc.sh'],
    entry_points={
        'console_scripts': [
            'TumOnc=TumOnc.__main__',
        ],
    },
    cmdclass={
        'publish': MinusPublish
    },
)
