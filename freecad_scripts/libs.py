#!/usr/bin/env python3
# Copyright (C) 2021, Miklos Maroti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys

# https://wiki.freecadweb.org/Embedding_FreeCAD
freecad_libs = [
    '/usr/local/lib/FreeCAD.so',
    '/usr/lib/freecad-python3/lib/FreeCAD.so',
    'c:/Program Files/FreeCAD 0.19/bin/FreeCAD.pyd',
]
for lib in freecad_libs:
    if os.path.exists(lib):
        path = os.path.dirname(lib)
        if path not in sys.path:
            sys.path.append(path)
        break
else:
    raise ValueError("FreeCAD library was not found!")

import FreeCAD                              # noqa
from FreeCAD import Units                   # noqa

femtools_libs = [
    '/usr/local/Mod/Fem/femtools',
    '/usr/share/freecad/Mod/Fem/femtools',
    'c:/Program Files/FreeCAD 0.19/Mod/Fem/femtools',
]
for lib in femtools_libs:
    if os.path.exists(lib):
        path = os.path.dirname(lib)
        if path not in sys.path:
            sys.path.append(path)
        path = os.path.abspath(os.path.join(lib, '..', '..'))
        if path not in sys.path:
            sys.path.append(path)
        path = os.path.abspath(os.path.join(lib, '..', '..', '..', 'Ext'))
        if path not in sys.path:
            sys.path.append(path)
        break
else:
    raise ValueError("femtools library was not found!")

from femtools.ccxtools import FemToolsCcx   # noqa
from femmesh.gmshtools import GmshTools     # noqa
