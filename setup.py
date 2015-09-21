import os
from setuptools import setup

setup(
    name="Pyabolism",
    version="0.4",
    author="Nick Fyson",
    author_email="mail@nickfyson.co.uk",
    description=("A module for constraint-based simulation of metabolic networks."),
    license = "BSD",
    keywords = "FBA metabolic network",
    # url = "http://pypi.python.org/pypi/pyabolism",
    packages=['pyabolism', 'pyabolism.simulate'],
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: BSD License",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2.7",
    ],
    test_suite='tests',
)
