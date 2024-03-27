from setuptools import setup

setup(
    name='pyfdstools',
    version='0.0.15',    
    description='This software is part of a python library to assist in developing and analyzing simulation results from Fire Dynamics Simulator (FDS). FDS is an open source computational fluid dynamics (CFD) software package developed by NIST. The sourcecode is available at: https://github.com/firemodels/fds',
    url='https://github.com/johodges/pyfdstools',
    author='Jonathan Hodges',
    author_email='johodges@vt.edu',
    license='MIT',
    packages=['pyfdstools'],
    
    include_package_data=True,
    install_requires=[  'libopencv>=0.0.1',
                        'matplotlib>=3.0',
                        'numpy>=1.17',
                        'numpy-stl>1.0',
                        'opencv-contrib-python>=4.7.0',
                        'pandas>=0.25',
                        'pyevtk>=1.0',
                        'scipy>=1.3.1',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
