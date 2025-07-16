from setuptools import setup, find_packages
import pathlib


# The directory containing this file
HERE = pathlib.Path(__file__).parent


# The text of the README file
README = (HERE / "README.md").read_text()


setup(name='nease',
      version='1.3.2',
      description='NEASE: A network-based approach for the enrichment of alternative splicing events',
      long_description=README,
      long_description_content_type="text/markdown",
      url='https://github.com/louadi/NEASE',
      author='Zakaria Louadi',
      author_email='Zakaria.louadi@tum.de',
      license='GPLv3',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3.6',
      install_requires=[
        'pandas>=1.3.0',
        'numpy>= 1.19',
        'networkx>=2.5',
        'plotly>=4.14.3',
        'scipy>=1.6.2',
        'statsmodels>=0.12.2',
        'gseapy>=0.10.7',
        'matplotlib>=3.4.1',
        'mygene>=3.2.2',
        'kaleido',
    ],

)
