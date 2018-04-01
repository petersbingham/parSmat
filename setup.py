# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'parSmat/README.md')

setup(name='parSmat',
      version='0.11',
      description='Python package to parametrise the multi-channel S-matrix using a pade approximation.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['parSmat'],
      package_data={'parSmat': ['tests/*', 'README.md']}
     )
