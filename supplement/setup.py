from setuptools import setup, find_packages

setup(
    name='dispersal_model',
    version='0.3',
    description='agent-based model of human dispersal',
    author='Gereon A. Kaipig',
    author_email='gereon.kaiping@geo.uzh.ch',
    url='https://github.com/Anaphory/dispersal-simulation',
    packages=['dispersal_model'],
    license="BSD (3 clause)",
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: BSD License',
    ],
    install_requires=[
        'cython',
        'attrs',
        'numpy',
        'matplotlib',
        'tifffile',
        'overpy',
        'shapely',
        'cartopy',
        'pyshp',
        'attrs',
        'h3',
    ],
    extras_require={
        'dev': [
            'mypy'
        ],
        'test': [
            'pytest>=3.6',
        ],
    },
    entry_points={
        'console_scripts': ['dispersal-model=dispersal_model.__main__:main'],
    },
    package_data={'dispersal_mdel': []},
)
