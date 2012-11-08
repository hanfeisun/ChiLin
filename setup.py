#!/usr/bin/env python

"""Description

"""

import sys
from setuptools import setup, find_packages

if sys.version < "2.6.0" or sys.version > "2.8.0":
    print "Please use a Python with higher version than 2.6.0"
    sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
    sys.exit(1)
    exit(1)

def main():
    setup(name="chilin",
          version="1.0.0",
          description="QC report pipline",
          author='Shenglin Mei, Qian Qin, Hanfei Sun',
          author_email='samleomei@gmail.com',
          packages = ["chilin"],
          package_dir={'chilin' : 'chilin/lib'},
          install_requires=['jinja2', 'argparse'],
          package_data = {'chilin':['db/*','template/*', 'conf/*', 'awk/*']},
          scripts=['chilin/scripts/ChiLin.py'],
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
