from distutils.core import setup


setup(name="xsgen",
      version="0.1",
      packages=["xsgen"],
      entry_points={"console_scripts": ["xsgen=xsgen.main:main"]})
