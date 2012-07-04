#!/usr/bin/env python
# Time-stamp: <2012-7-07 17:24:20 Shenglin Mei>

"""Description

"""

import os
import sys
#from setuptools import setup, find_packages
from distutils.core import setup, Extension
try:
    import py2exe
except ImportError:
    pass
try:
    import py2app
except ImportError:
    pass

def setup_conf():
	path=os.getcwd()
	


def main():
	path=os.getcwd()
	print path
	if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
		sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
		sys.exit(1)
	setup(name="QCreport",
		version="1.0.0",
		description="QC report pipline",
		author='Shenglin Mei',
		author_email='samleomei@gmail.com',
		package_dir={'QCreport' : 'lib'},
		packages=['QCreport'],
		scripts=['bin/QCreport'],
		console=['bin/QCreport'],
		app=['bin/QCreport'],
		classifiers=[
		'Development Status :: 5 - productive',
		'Environment :: Console',
		'Intended Audience :: Developers',
		'License :: OSI Approved :: Artistic License',
		'Operating System :: MacOS :: MacOS X',
		'Operating System :: Microsoft :: Windows',
		'Operating System :: POSIX',
		'Programming Language :: Python',
		],
		)
if __name__ == '__main__':
	main()
