#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import re
from setuptools import setup


# load version form _version.py
VERSIONFILE = "spe_parser/_version.py"
with open(VERSIONFILE, 'rt', encoding='utf-8') as f:
    verstrline = f.read()

VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError(f"Unable to find version string in {VERSIONFILE}.")

# README
with open('README.md', 'r', encoding='utf-8') as f:
    long_desc = f.read()

# module

setup(name='spe_parser',
      version=verstr,
      author="Keisuke Fujii",
      author_email="fujiisoup@gmail.com",
      description=("Python package to read spe file."),
      long_description=long_desc,
      long_description_content_type='text/markdown',
      license="BSD 3-clause",
      keywords="imaging, Princeton Instruments",
      url="https://github.com/fujiisoup/spe_parser",
      packages=["spe_parser"],
      package_dir={'spe_parser': 'spe_parser'},
      py_modules=['spe_parser.__init__'],
      test_suite='testings',
      install_requires=[
        'numpy>=1.10',
        'untangle'],
      classifiers=['License :: OSI Approved :: BSD License',
                   'Natural Language :: English',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Microsoft :: Windows',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 3.6',
                   'Topic :: Scientific/Engineering :: Physics'],
      )