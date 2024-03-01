from setuptools import setup, find_packages

setup(
    name='SigConfide',
    version='0.1.0',
    author='Marcin Wierzbiński, Damian Wójtowicz',
    author_email='m.wierzbinski@uw.edu.pl',
    description='Analysis and decomposition of cancer signatures using bootstrap and cross-validation methods',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/marcin119a/SigConfide',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'quadprog'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        #'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)