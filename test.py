#!/usr/bin/env python

import os, sys
import nose

if not nose.run(argv=[__file__, '-v','--with-doctest']): sys.exit()

isolated = ("nsfd2o", )

for p in isolated:
    ret = os.system(" ".join( (sys.executable, "test/test_%s.py"%p) ))
    if ret != 0: sys.exit()

ret = os.system("cd doc/sphinx && make doctest")
if ret != 0: sys.exit()