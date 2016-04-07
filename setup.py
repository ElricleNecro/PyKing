#! /usr/bin/env python
# -*- coding:Utf8 -*-

import os
import numpy as np
try:
	import commands
except:
	import subprocess as commands

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def pkgconfig(*packages, **kw):
	flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
	for token in commands.getoutput("pkg-config --libs --cflags %s" % ' '.join(packages)).split():
		kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
	return kw


src_files = [ "cKing.pyx" ]

opt = pkgconfig("king")
opt["include_dirs"] += [
		'.',
		np.get_include(),
	] #'/home/plum/.local/lib/python3.3/site-packages/Cython/Includes/']

print(opt)

setup(
	name         = 'King',
	data_files   = [(os.path.join(opt["include_dirs"][0], 'king'), ['King.pxd'])],
	version      = '1.0',
	description  = 'Wrapper cython around my king library.',
	author       = 'Guillaume Plum',
	cmdclass     = {'build_ext': build_ext},
	ext_modules  = [
		Extension("King", src_files, **opt)
	]
)

