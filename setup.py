# -*- coding: utf-8 -*-

from distutils.core import setup
import os
import shutil
shutil.copy('README.md', 'parsmat/README.md')

dir_setup = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(dir_setup, 'parsmat', 'release.py')) as f:
    # Defines __version__
    exec(f.read())

setup(name='parsmat',
      version=__version__,
      description='Python package to parametrise the multi-channel S-matrix using a pade approximation.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['parsmat'],
      package_data={'parsmat': ['tests/*', 'README.md']}
     )
