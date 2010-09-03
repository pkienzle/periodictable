#!/usr/bin/env python

import sys
import os
for f in ['core','covalent_radius','crystal_structure',
          'density','formulas','magnetic_ff','mass',
          'nsf','private','xsf']:
    os.system(sys.executable+" "+f+".py")
