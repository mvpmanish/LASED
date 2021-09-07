from setuptools import setup, find_packages

VERSION = '0.3' 
DESCRIPTION = 'LASED'
LONG_DESCRIPTION = 'A Laser-Atom Interaction Simulator using Quantum Electrodynamics'

with open("../README.md", "r") as fh:
    long_description = fh.read()

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
            "sympy >= 1.8"
        ], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'laser-atom', 'simulation', 'quantum', 'quantum electrodynamics', 'physics'],
        classifiers= [
            "Development Status :: 1 - Planning",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
            "Natural Language :: English"
        ]
)