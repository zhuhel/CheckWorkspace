from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from setuptools import find_packages
from setuptools import setup

description="A light script to check workspaces"

setup(
    name="checkws",
    version="0.0.1",
    description=description,
    author="Xiangyang Ju",
    packages=find_packages(),
    url='https://gitlab.cern.ch/xju/CheckWorkspace',
    install_requires=[
        "argparse",
        "numpy"
    ],
    scripts=[
        'scripts/check_ws',
    ]
)
