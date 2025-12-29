import glob
import os
import re
from setuptools import setup, find_packages

# Read version from __init__.py
def get_version():
    init_path = os.path.join('pangenomium', '__init__.py')
    with open(init_path) as f:
        for line in f:
            if line.startswith('__version__'):
                return re.search(r'"([^"]+)"', line).group(1)
    raise RuntimeError("Unable to find version string")

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pangenomium',
    version=get_version(),
    description='Dereplicate genomes and proteins into pangenomes and orthologs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Josh L. Espinoza',
    author_email='jespinoz@jcvi.org',
    url='https://github.com/jolespin/pangenomium',
    license='Apache-2.0',
    
    # Automatically find packages
    packages=find_packages(),
    
    # Dependencies
    python_requires='>=3.6',
    install_requires=requirements,
    
    # Entry points - creates the 'pangenomium' command
    entry_points={
        'console_scripts': [
            'pangenomium=pangenomium.cli:main',
        ],
    },
    
    # Additional scripts
    scripts=glob.glob('scripts/*.py'),
    
    # Include non-Python files
    include_package_data=True,
    
    # Classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    
    keywords='bioinformatics genomics pangenome clustering',
)
