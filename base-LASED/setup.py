from setuptools import setup, find_packages
from pybind11 import get_cmake_dir # Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext

VERSION = '1.0.1'
DESCRIPTION = 'LASED'
LONG_DESCRIPTION = 'A Laser-Atom Interaction Simulator using Quantum Electrodynamics'

with open("../README.md", "r") as fh:
    long_description = fh.read()

# Setup the C++ module
ext_modules = [
    Pybind11Extension("LASED.CppLASED",
        ["src/main.cpp"],
        include_dirs = ["include"],
        define_macros = [('VERSION_INFO', VERSION)],
        ),
]

# Setting up
setup(
       # the name must match the folder name 'LASED'
        name = "LASED",
        version = VERSION,
        author = "Manish Patel",
        author_email="<mvpmanish@gmail.com>",
        description = DESCRIPTION,
        long_description = long_description,
        long_description_content_type = "text/markdown",
        packages=find_packages(),
        install_requires=[
            "numpy >= 1.20",
            "scipy >= 1.6.0",
            "sympy >= 1.8",
            "pybind11 >= 2.8.1" # package that resuqired for the cpp module
        ], # add any additional packages that
        # needs to be installed along with your package. Eg: 'caer'
        ext_modules=ext_modules, # Setup the C++ module
        keywords=['python', 'laser-atom', 'simulation', 'quantum', 'quantum electrodynamics', 'physics'],
        classifiers= [
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
            "Natural Language :: English",
            "Topic :: Scientific/Engineering :: Physics"
        ]
)
