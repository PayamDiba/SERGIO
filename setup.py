import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '1.9.2'
PACKAGE_NAME = 'SERGIO-scSim'
AUTHOR = 'Payam Dibaeinia'
AUTHOR_EMAIL = 'dibaein2@illinois.edu'
URL = 'https://github.com/PayamDiba/SERGIO'

LICENSE = 'GNU GENERAL PUBLIC LICENSE'
DESCRIPTION = 'SERGIO is a simulator for single-cell expression data guided by gene regulatory networks.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'pandas',
      'absl-py',
      'networkx',
      'cma',
      'matplotlib',
      'scikit-learn'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      include_package_data=True,
      packages=find_packages('SERGIO'),
      python_requires='>3.5.2',
      classifiers =(
            "Programming Language :: Python :: 3", ),
      )
