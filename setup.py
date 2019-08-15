"""
@author   Thomas B. Wild
@email:   twild@umd.edu or tombernardwild@gmail.com
@Project: pysedsim 1.0.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files
"""


class VersionError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Must have setuptools installed to run setup.py. Please install and try again.")
    raise


def readme():
    with open('README.md', encoding='utf-8') as f:
        return f.read()


def get_requirements():
    with open('requirements.txt') as f:
        return f.read().split()


setup(
    name='pysedsim',
    version='1.0.0',
    packages=find_packages(),
    url='https://github.com/FeralFlows/pysedsim',
    license='BSD 2-Clause',
    author='Thomas B. Wild; Abigail N. Birnbaum; Patrick M. Reed; Daniel P. Loucks',
    author_email='twild@umd.edu',
    description='An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting, Design, and Operation Alternatives',
    long_description=readme(),
    install_requires=get_requirements(),
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    classifiers=[
        "Programming Language :: Python :: 2",
        "Operating System :: OS Independent"]
)
