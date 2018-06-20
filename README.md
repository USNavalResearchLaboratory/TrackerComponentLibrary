**Tracker Component Library Release 3.5, May 2018**
https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary

A paper describing a number of features of the library is<br>
D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid
Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, no.
5, pp. 18-27, May. 2017.

These are the release notes for the version 3.5 of the Tracker Component
Library. The Tracker Component Library is a collection of Matlab routines
(some written in C/C++ requiring compilation compilation) for simulating
and tracking targets in various scenarios. Due to the complexity of the
target tracking problem, a great many routines can find use in other areas
including combinatorics, astronomy, and statistics.

As of version 3.0, the library has been split into two parts. This is the
publicly available part. A number of functions have been placed into a
supplement that is only available to the U.S. Government Agencies and their
contractors. More advanced algorithms added to the library will be placed
in the limited distribution supplement. The limited distribution supplement
can be downloaded by authorized individuals from di2e.net at:
https://bitbucket.di2e.net/projects/TCL

Those looking for magnetic field sythesis code might want to look at
./Sample Code/Magnetic Models/ .

Those looking to get a quick idea of a very simple end-to-end tracking
algorithm (track initiation, data association, maintenance, and
termination) might want to look at
./Sample Code/Basic Tracking Examples/demo2DIntegratedDataAssociation.m
to see an example of a complete GNN-JIPDAF tracker run on a simple 2D
tracking scenario, which can be easily transformed into a normal JIPDAF
by changing a single parameter. The tracker is made from the modular
components of the library and adapting the tracker to 3D and other
scenarios is as simple as changing the proper components. A lot of other
sample code is also provided to demonstrate the use of other parts of the
library.

To use the full library, add the library and all of its subfolders to your
active path in Matlab. Some functions are available as C/C++ files for use
in Matlab either because they use third-party libraries (and must be
compiled to be used) or because the native Matlab implementations provided
are too slow in certain circumstances. These files can be compiled by
running the CompileCLibraries function. Precompiled code is not distributed
with the library. Note that a C/C++ compiler supported by Matlab must be
installed. See below for comments regarding compilation.

NOTABLE CHANGES SINCE VERSION 3.01

- More compiled versions of 2D and 3D assignment code are available and the
  C implementations have been more modularized to make then simpler to
  remove from the Matlab parts of the library and use elsewhere.
- Sped up the dense interior point linear programming method by a large
  factor compared to the textbook implementation.
- A few new optimization methods.
- Improved the probability distributions. Added new ones, such as the
  wrapped normal distribution, and added entropy to many of them.
- For the many functions that have no Matlab implementation and must be
  compiled to be used (mostly those using third-party libraries), shell .m
  files have been added to the library, so that one can type "help" and the
  name of the function to get information on how to use it. Also, if the
  function is called without having compiled it, then an error alerts the
  user to the need to compiled the function.
- Bug fixes and small changes.
- Updated third-party code. Removed reliance on third-party functions for
  line searches.

NOTABLE CHANGES SINCE VERSION 3.0

- Minor corrections.
- Minor modifications to files to support the final release of the limited
  distribution supplement to the library.
- New 2D and 3D assignment algorithms, among other additions.

NOTABLE CHANGES SINCE VERSION 2.5

- Many new functions have been added.
- Many corrections have been made and a number of functions in different
  areas have been added or made more versatile/ stable.
- A number of routines related to signal processing have been added to the
  library in the Mathematical Functions/Signal Processing folder.
- Less reliance on third-party libraries and compiled code. NASA's SPICE
  library has been removed and the function readJPLEphem can be used to
  read in NASA JPL's ephemerides without compiling anything. Additionally,
  the directGeodeticProb function has been rewritten in Matlab and does not
  use GeographicLib, though the GeographicLib library is still used to
  solve the indirect geodetic problem.
- Spherical harmonic synthesis has been simplified and now works with
  vector coefficients and complex coefficients. Additionally, in addition
  to obtaining a potential and a gradient, one can obtain a Hessian matrix,
  opening the door to simulation of gravity gradiometer measurements.
  Fitting of low to moderate order spherical harmonic coefficients to
  measured data is now possible with the fitSpherHarmonics function.
- Arbitrary-order multivariate B-spline fitting, evaluation,
  differentiation, and integration has been added. 
- Linear dynamic models expressed in non-Cartesian coordinate systems have
  been added along with conversions of states between coordinate systems.
  See /Dynamic Models/Continuous Time/Non-Cartesian Position/ and
  /Coordinate Systems/State Coordinate System Conversion/

COMPILED CODE:

The compilation of the library has been tested under Matlab2018a for Mac OS
X 10.11 using the Matlab-supported compiler offered free by Apple. Under
Windows 10, compilation of the library was verified using Microsoft Visual
Studio Professional 2017 and also using minGW 5.3 as the compiler.
Precompiled code is not distributed with the library.

EXTERNAL SOLVERS:

The folder
./Mathematical Functions/Polynomials/Generic Multivariate Polynomials/
contains the function solvePolySysWithExtProg. This is an interface to
external solvers for simultaneous multivariate polynomial systems. The
solvers are external programs and are not included. Though the functions
polyRootsMultiDim and solveQuadBivarEq can solver certain types of
simultaneous multivariate polynomials, they are limited. For more general
systems, an external solver is the best choice. The functions
polyMeasConvert, polyMeasConvertAsync, and DopplerOnlyInit6D are located in
subfolders entitled 'Uses External Solver" and use the function
solvePolySysWithExtProg. Supported external solvers for the function
solvePolySysWithExtProg are PHCpack, Bertini, and the certified solver in
Macaulay2. Note that Macaulay2 tends to fail more often than the other
solvers.

The solvers can be downloaded from:
PHCpack for Mac OS X, Windows and Linux:
http://homepages.math.uic.edu/~jan/download.html
Bertini (version 1.5) for Mac OS X and Linux and Windows, only with Cygwin:
https://bertini.nd.edu/download.html
Macaulay2 for Mac OS X, Linux, and Windows:
http://www.math.uiuc.edu/Macaulay2/

For the external solvers to work, they need to be added to the search path
used in Matlab so that the phs/bertini/m2 commands can be called using the
system command in Matlab.

DATA SOURCES:

If the library was downloaded without associated data files for solar
system ephemerides as well as magnetic, gravitational and terrain models,
then the relevant data files will have to be downloaded The original 
sources of the data are:

1) Earth2014 terrain model:
https://geodesy.curtin.edu.au/research/harmonic/
The files Earth2014.BED2014.degree2160.bshc,
Earth2014.ICE2014.degree2160.bshc, Earth2014.RET2014.degree2160.bshc,
Earth2014.SUR2014.degree2160.bshc, Earth2014.TBI2014.degree2160.bshc should
be placed into ./Terrain/data . Though data files to degree and order
10,800 are also available, they are not currently supported as the function
spherHarmonicEval does not use extended precision arithmetic, which is
required to avoid overflows when using the model.

2) EGM2008 terrain model (The DTM2006.0 model):
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html
The file Coeff_Height_and_Depth_to2190_DTM2006.0 should be placed into 
./Terrain/data . The first time using the getEGM2008TerrainCoeffs function
will be significantly slower than subsequent calls, assuming the full model
is loaded the first time.

3) The EMM2017 magnetic field coefficients:
https://www.ngdc.noaa.gov/geomag/EMM/
The coefficients are included with the software that is provided. The files
having names like EMM2015.COF and EMM2000.COF should all be zipped into an
archive with the name EMM_Coefficients.zip and placed in ./Magnetism/data .
The files with names like EMM2015SV.COF andEMM2000SV.COF should all be
zipped into an archive with the name EMM_Secular_Variations.zip and placed
into ./Magnetism/data . If all of the data for a single year is loaded,
then a .mat file containing the data is placed in  the data folder to speed
up subsequent calls to the function.

4) The IGRF12 magnetic field coefficients:
http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
The file http://www.ngdc.noaa.gov/IAGA/vmod/igrf12coeffs.txt can be
downloaded and placed into the folder ./Magnetism/data .

5) The WMM magnetic field coefficients:
https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml
The coefficient file WMM2015COF.zip shown be downloaded, unzipped and the
file WMM.COF placed into ./Magnetism/data .

6) EGM2008 gravitational coefficients
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
Download the file  EGM2008_to2190_TideFree.gz, which is under "EGM2008 Tide
Free Spherical Harmonic Coefficients". Decompress the file and place it in
./Gravity/data .

7) EGM96 gravitational coefficients
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
The file egm96.z should be downloaded and decompressed, resulting in a file
names egm96 without an extension, which should be placed in
./Gravity/data .

8) EGM2008 parameters needed to convert height anomalies to geoid
   undulations:
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
The "Correction Model" file Zeta-to-N_to2160_egm2008.gz should be
downloaded, ungzipped and placed in ./Gravity/data .

9) EGM96 parameters needed to convert height anomalies to geoid
   undulations:
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
The file corrcoef.z should be decompressed and placed in ./Gravity/data .

10) NASA JPL's DE430t planetary/solar/lunar ephemerides
ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t
The file linux_p1550p2650.430t should be downloaded and placed in 
./Astronomical Code/data
for use by the readJPLEphem function. The readJPLEphem function can read
other NASA ephemerides, but the DE430t set are the default.

11) NASA JPL's GL0900C gravitational model of the Moon
http://pds-geosciences.wustl.edu/missions/grail/default.htm
Download the file jggrx_0900c_sha.tab and place it in ./Gravity/data .

12) The FES2004 tide model of the effects of ocean Earth tides.
http://maia.usno.navy.mil/conv2010/convupdt/convupdt_c6.html
The file fes2004_Cnm-Snm.dat should be downloaded and placed in
./Gravity/Tides/data .

14) The Hipparcos 2 Star Catalog
http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/311
Place the file hip2.dat in ./Astronomical Code/data .

15) A low-resolution image of the "Blue Marble" Earth from NASA with land,
    ocean color, sea ice and clouds.
https://visibleearth.nasa.gov/view.php?id=57735
Place the file in Misc/data . This is used as the default map to show on
the Earth given by the plotMapOnEllipsoid function.

May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.<br>
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

LICENSE:

The source code is in the public domain and not licensed or under
copyright. The information and software may be used freely by the public.
As required by 17 U.S.C. 403, third parties producing copyrighted works
consisting predominantly of the material produced by U.S. government
agencies must provide notice with such work(s) identifying the U.S.
Government material incorporated and stating that such material is not
subject to copyright protection.

Derived works shall not identify themselves in a manner that implies an
endorsement by or an affiliation with the Naval Research Laboratory.

RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
OF RECIPIENT IN THE USE OF THE SOFTWARE.