from setuptools import setup, find_packages

VERSION = '0.1' 
DESCRIPTION = 'LASED'
LONG_DESCRIPTION = 'A Laser-Atom Interaction Simulator using Quantum Electrodynamics'

# Setting up
setup(
       # the name must match the folder name 'LASED'
        name = "LASED", 
        version = VERSION,
        author = "Manish Patel",
        author_email="<mvpmanish@gmail.com>",
        description = DESCRIPTION,
        long_description = LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'laser-atom', 'simulation', 'quantum', 'quantum electrodynamics', 'physics'],
        classifiers= [
            "Development Status :: 1",
            "Intended Audience :: Education/Research",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: Linux :: Ubuntu",
            "Operating System :: Microsoft :: Windows",
        ]
)