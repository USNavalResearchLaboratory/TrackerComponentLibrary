External libraries that are called by functions within the Tracker
Component Library are kept here. In many instances, only a few of the
routines in the library are used, though with the exception of the SGP4
library, the entire libraries are provided here. None of the parts of the
external libraries that are used is subject to copyleft restrictions,
though they have varying licenses. The Splitting Conic Solver Library and
the Eigen library can use some routines that rely are under an LGPL
(copylefted) license. Those routines are not used and have been placed in a
zip archive in the respective folders. The Eigen library is additionally
compiled with an option that prevents the use of its LGPLd code.

The sources of the external libraries are as follows:

1) The dividing rectangles algorithm for global optimization is taken from
   the NLOpt library at
https://http://ab-initio.mit.edu/wiki/index.php/NLopt

2) LibAIS for decoding data from ships' transponders:
https://github.com/schwehr/libais

3) Liblbfgs for optimization
http://www.chokkan.org/software/liblbfgs/

4) The C Implementation of the NRLMSISE-00 empirical atmospheric model:
http://www.brodo.de/space/nrlmsise/

5) The Splitting Conic Solver:
The Matlab version (used here ) is from
https://github.com/bodono/scs-matlab
but a C version is also available from:
https://github.com/cvxgrp/scs

6) The implementation of the SGP4 orbital propagator:
http://www.centerforspace.com/downloads/

7) The International Astronomical Union's Standards of Fundamental
   Astronomy Library:
http://www.iausofa.org/index.html

8) The Eigen Library (Version 3.4.0 is used):
https://eigen.tuxfamily.org/
