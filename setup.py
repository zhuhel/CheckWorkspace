from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from setuptools import find_packages
from setuptools import setup

description="Script to check workspaces"

setup(
    name="checkws",
    version="0.1.0",
    description=description,
    author="Xiangyang Ju",
    packages=find_packages(),
    url='https://gitlab.cern.ch/xju/CheckWorkspace',
    install_requires=[
        "argparse",
        "numpy",
        'root_plot_utils',
    ],
    dependency_links=[
        'https://github.com/xju2/root_plot_utils/tarball/v1.0.1#egg=root_plot_utils-1.0.1'
        # 'https://github.com/xju2/root_plot_utils.git'
    ],
    # setup_requires=['root_plot_utils==1.0.1'],
    scripts=[
        'scripts/check_ws',
        'scripts/print_np',
        'scripts/draw_pulls',
    ]
)
