#!/usr/bin/env python

from distutils.core import setup

setup(name='phylolabel',
      version='0.1',
      description='Label higher-order taxa in a phylogeny using a reference taxonomy',
      author='Ben Morris',
      author_email='ben@bendmorris.com',
      url='https://github.com/bendmorris/phylolabel',
      packages=['phylolabel'],
      package_dir={
                'phylolabel':''
                },
      entry_points={
        'console_scripts': [
            'treestore = phylolabel.phylolabel:main',
        ],
      },
      )