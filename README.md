**Tracker Component Library Release 3.0, September 2017**
https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary

A paper describing a number of features of the library is<br>
D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid
Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, no.
5, pp. 18-27, May. 2017.

These are the release notes for the version 3.0 of the Tracker Component
Library. The Tracker Component Library is a collection of Matlab routines
(some written in C/C++ requiring compilation compilation) for simulating
and tracking targets in various scenarios. Due to the complexity of the
target tracking problem, a great many routines can find use in other areas
including combinatorics, astronomy, and statistics.

As of version 3.0, the library has been split into two parts. This is the
publicly available part. A small number of functions have been placed
into a supplement that is only available to the U.S. DoD and DoD and U.S.
DoD contractors (Distribution Statement D), though certain parts may be
available to other U.S. Government Agencies (DIstributions Statement C)
upon request. Many future additions to the library will be placed in the
the limited distribution supplement. Instructions for obtaining the limited
distribution potion of the library will be posted in this readme after the
supplement has been made available.

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

The compilation of the library has been tested under Matlab2017a for Mac OS
X 10.11 using clang-800.0.42.1. However, the code should work under other
versions of Matlab, compilers, and operating systems.  To use the full
library, add the library and all of its subfolders to your active path in
Matlab. Some functions are available as C/C++ files for use in Matlab
either because they use third-party libraries or because the native Matlab
implementations provided are too slow to use. If the library was downloaded
without precompiled versions of these files, then they can be compiled by
running the CompileCLibraries function. Note that a C/C++ compiler
supported by Matlab must be installed. Additionally, to decode AIS data,
the compiler must support C++11. As the commands for the compilation of the
AIS data are at the end of the CompileCLibraries function, if a compiler
without C++11 support is used, the rest of the library should be
successfully compiled before an error occurs while compiling the AIS code.

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

September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.