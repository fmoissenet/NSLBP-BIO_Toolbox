# -*- coding: utf-8 -*-
from setuptools import setup,find_packages
import os,sys
import logging

developMode = False
if len(sys.argv) >= 2:
    if sys.argv[1] == "develop": developMode = True
if developMode:
    logging.warning("You have sleected a developer model ( local install)")


VERSION ="1.0.0"


#------------------------- INSTALL--------------------------------------------
setup(name = 'pyCGM2-flow',
    version = VERSION,
    author = 'Fabien Leboeuf',
    author_email = 'fabien.leboeuf@gmail.com',
    description = "pyCGM2 extension for working with configurable flow",
    long_description= "",
    url = 'https://pycgm2.github.io',
    keywords = 'python CGM CGM2 Vicon PluginGait',
    packages=find_packages(),
	include_package_data=True,
    license='CC-BY-SA',
	install_requires = ['numpy<1.17.0',
                        'scipy==1.2.1',
                        'matplotlib<3.0.0',
                        'pandas ==0.19.1',
                        'enum34>=1.1.2',
                        'configparser>=3.5.0',
                        'beautifulsoup4>=3.5.0',
                        'pyyaml>=3.13.0',
                        'yamlordereddictloader>=0.4.0',
                        'xlrd >=0.9.0',
                        'lxml>=4.4.1',
                        'openpyxl>=2.6.3',
                        'xlwt>=1.3.0',
                        'pytest==4.6.5',
                        'pymongo==3.11.0',
                        'Jinja2==2.11.2'],
    #'qtmWebGaitReport>=0.0.1'],
    classifiers=['Programming Language :: Python',
                 'Programming Language :: Python :: 2.7',
                 'Operating System :: Microsoft :: Windows',
                 'Natural Language :: English'],
    #scripts=gen_data_files_forScripts("Apps/ViconApps")
    entry_points={
          'console_scripts': [
                # NEXUS
                'pyCGM2f-init  =  pyCGM2f.commands:initiatingCommand',
                'pyCGM2f-edit  =  pyCGM2f.commands:editingCommand',
                'pyCGM2f-process  =  pyCGM2f.commands:processingCommand',
                'pyCGM2f-report  =  pyCGM2f.commands:reportingCommand',
                'pyCGM2f-check  =  pyCGM2f.commands:checkingCommand',
                'pyCGM2f-push  =  pyCGM2f.commands:pushingCommand'
          ]
      },
    )
