import setuptools
from distutils.core import setup
from pathlib import Path

REQUIREMENTS_FILE = 'requirements.txt'
requirements = Path(REQUIREMENTS_FILE).read_text().splitlines()

setup(
    name='Needal',
    version='1.0',
    description='Annotation tool for Text-Fabric',
    author='Cody Kingham',
    license='MIT',
    url='https://github.com/codykingham/needal',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(include=['needal']),
    install_requires=requirements,
)

