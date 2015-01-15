#! /usr/bin/env python
try:
    from setuptools import setup
    has_setuptools = True
except ImportError:
    has_setuptools = False
    from distutils.core import setup

if has_setuptools:
    setup_kwargs = {
        "entry_points": {"console_scripts": ["xsgen=xsgen.main:main"]}
    }
else:
    setup_kwargs = {
        "scripts": ["scripts/xsgen"]
    }

setup(name="xsgen",
      version="0.1",
      packages=["xsgen"],
      **setup_kwargs)
