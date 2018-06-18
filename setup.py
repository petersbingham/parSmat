# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'parsmat/README.md')

setup(name='parsmat',
      version='1.0.0',
      description='Python package to parametrise the multi-channel S-matrix using a pade approximation.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['parsmat'],
      package_data={'parsmat': ['tests/*', 'README.md']}
     )
