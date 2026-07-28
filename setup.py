"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-chem-netmhcpan',
    version='0.2.0',
    description='Scipion framework plugin for MHC-I/MHC-II promiscuity prediction with NetMHCpan-4.2 and NetMHCIIpan-4.3',
    long_description=long_description,
    url='https://github.com/Lvera-code/scipion-chem-netmhcpan',
    author='Enzo Sierra',
    author_email='enzogael57@gmail.com',
    keywords='scipion epitope mhc netmhcpan netmhciipan cytotoxic t-cell t-helper',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    entry_points={
        'pyworkflow.plugin': 'netmhcpan = netmhcpan'
    }
)
