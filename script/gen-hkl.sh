#!/bin/bash

../bin/gen-hkl XYZIN 1000006.pdb HKLOUT 1000006.mtz \
<<eof-gen-hkl
TITL generate structure factors
RESO 1.0
BFAC 5.0
EGAU 1.2
END
eof-gen-hkl
