
from setuptools import setup, find_packages

setup(
    name='imgCIF_Creator',
    version='0.0.1',
    packages=find_packages(include=['imgCIF_Creator', 'imgCIF_Creator.*']),
    install_requires=[
        'PyYAML',
        'numpy',
        'click',
        'PyCifRW',
        'h5py',
        'cbf',
    ],
    entry_points={
        "console_scripts": [
            "creator.py=imgCIF_Creator.creator:main",
        ],
    },
    include_package_data=True,

)
