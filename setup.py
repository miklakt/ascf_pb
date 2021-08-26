import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

setup(
    name='ascf_pb',
    version='0.0.5',
    description='Analytical Self-Consistent Field Polymer Brushes',
    author='Laktionov Mikhail',
    author_email = 'miklakt@gmail.com',
    packages=['ascf_pb', 'ascf_pb.topology'],
    install_requires=['numpy', 'scipy'],
    #entry_points={
    #    "console_scripts": [
    #        "ascf_pb=ascf_pb.__main__:main",
    #    ]
    #},
)