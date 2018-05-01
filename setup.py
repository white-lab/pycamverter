
import os
from setuptools import setup, find_packages

from pycamv import __version__


REQUIREMENTS_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), "requirements.txt",
    )
)


setup(
    name="pycamverter",
    version=__version__,
    description=(
        "Utility for converting searched mass spec data into a format "
        "readable by CAMV"
    ),
    url="https://github.com/white-lab/pycamverter",
    author="Nader Morshed",
    author_email="morshed@mit.edu",
    license="BSD",
    packages=find_packages(exclude=["*.tests", "tests"]),
    install_requires=[
        "numpydoc>=0.7",
        "openpyxl>=2.5.0",
        "pymzML",
        "requests>=2.18.4",
        "scipy>=1.0.1",
    ],
    dependency_links=[
        "git+git://github.com/naderm/pymzML.git",
    ],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        # "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering",
    ],
    test_suite="tests",
)
