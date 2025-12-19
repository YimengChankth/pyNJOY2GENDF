# Tutorial
# https://github.com/areed1192/sigma-coding for help with setup
# https://www.youtube.com/watch?v=avtVgeO7DV8&list=PLcFcktZ0wnNm_P4SDR4aDVGqt4MeuSMfS&index=3 
# Install in editable mode, so any changes in this folder will be made immediately instead of using uninstall and install everytime changes are made 
# pip3 install -e .


from setuptools import setup
from setuptools import find_namespace_packages

setup(
    # Define library name, this is what is used along with pip install
    name='pyNJOY2GENDF',

    # Define the author of the repository
    author='Yimeng Chan',

    # Define the author's email
    author_email='yi.chan@psi.ch',

    # Link to url
    # url=

    # Define the version of this library
    # Read as:
    # Major version: 0
    # Minor version: 0
    # Maintanence version: 0
    # .dev1: development release
    # a1: alpha release
    # b1: beta release

    version='0.0.0a1',

    description='A python based wrapper to write NJOY files to produce multi-group ENDF (GENDF) files.',

   # not strictly for this python version
    python_requires='>=3.13.9', 

    # packages I want to "build"
    packages=['pyNJOY2GENDF'],
    
    install_requires=[
        "numpy>=2.3.3", # not strictly necessary but untested on ver 1.
        "tabulate",
    ],

)

