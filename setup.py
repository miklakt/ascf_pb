import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

setup(
    name='ascf_pb',
    version='0.0.2',
    description='Analytical Self-Consistent Field for planar Polymer Brush',
    author='Laktionov Mikhail',
    author_email = 'miklakt@gmail.com',
    packages=['ascf_pb'],
    install_requires=['numpy', 'scipy'],
    entry_points={
        "console_scripts": [
            "ascf_pb=ascf_pb.__main__:main",
        ]
    },
)