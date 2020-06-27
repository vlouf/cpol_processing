#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = "cpol_processing"
DESCRIPTION = """Radar PPIs data processing, quality control, filtering, attenuation
correction, dealiasing, unfolding, hydrometeors calculation, rainfall rate estimation."""
URL = "https://github.com/vlouf/cpol_processing"
EMAIL = "valentin.louf@bom.gov.au"
AUTHOR = "Valentin Louf"

# What packages are required for this module to be executed?
REQUIRED = [
    "arm_pyart",
    "numpy",
    "csu_radartools",
    "crayons",
    "netCDF4",
    "scipy",
    "unravel",
]

here = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = "\n" + f.read()

class PublishCommand(Command):
    """Support setup.py publish."""

    description = "Build and publish the package."
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print("\033[1m{0}\033[0m".format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status("Removing previous builds…")
            rmtree(os.path.join(here, "dist"))
        except FileNotFoundError:
            pass

        self.status("Building Source and Wheel (universal) distribution…")
        os.system("{0} setup.py sdist bdist_wheel --universal".format(sys.executable))

        self.status("Uploading the package to PyPi via Twine…")
        os.system("twine upload dist/*")

        sys.exit()


setup(
    name=NAME,
    version='2.6.2',
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={"cpol_processing": ["data/GM_model_CPOL.pkl.gz"]},
    install_requires=REQUIRED,
    include_package_data=True,
    license="ISC",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    keywords="radar weather meteorology dual-polarization hydrometeors rainfall",
    cmdclass={"publish": PublishCommand,},
)
