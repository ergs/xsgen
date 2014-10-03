try:
    from setuptools import setup
    has_setuptools = True
except ImportError:
    from distutils.core import setup

if has_setuptools:
    setup_kwargs = {
        "entry_points": {"console_scripts": ["xsgen=xsgen.xsgen:main"]}
    }
else:
    setup_kwargs = {
        "scripts": ["xsgen/xsgen.py"]
    }

setup(name="xsgen",
      version="0.1",
      packages=["xsgen"],
      **setup_kwargs)
