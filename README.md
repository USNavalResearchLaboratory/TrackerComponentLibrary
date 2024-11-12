**Tracker Component Library Release 6.1, September 2024**
https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary

A paper describing a number of features of the library is<br>
D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid
Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, no.
5, pp. 18-27, May. 2017.

Those looking for magnetic field synthesis code might want to look at
./Sample Code/Magnetic Models/ .

These are the release notes for the version 6.1 of the Tracker Component
Library. The Tracker Component Library is a collection of Matlab routines
for simulating and tracking targets in various scenarios. Due to the
complexity of the target tracking problem, a great many routines can find
use in other areas including combinatorics, astronomy, and statistics.

Making this code available and open source is in the spirit of the Federal
Source Code Policy, which is OMB M-16-21, which is given in
https://www.whitehouse.gov/wp-content/uploads/legacy_drupal_files/omb/memoranda/2016/m_16_21.pdf 
https://code.gov/

Those looking to get a quick idea of a very simple end-to-end tracking
algorithm (track initiation, data association, maintenance, and
termination) might want to look at
./Sample Code/Basic Tracking Examples/demo2DIntegratedDataAssociation.m
to see an example of a complete GNN-JIPDAF tracker run on a simple 2D
tracking scenario, which can be easily transformed into a normal JIPDAF
by changing a single parameter. The tracker is made from the modular
components of the library. A lot of other sample code is also provided to
demonstrate the use of other parts of the library.

To use the full library, add the library and all of its subfolders to your
active path in Matlab. Some functions are available as C/C++ files for use
in Matlab either because they use third-party libraries (and must be
compiled to be used) or because the native Matlab implementations provided
are too slow in certain circumstances. These files can be compiled by
running the CompileCLibraries function. Precompiled code is not distributed
with the library. Note that a C/C++ compiler supported by Matlab must be
installed. See below for comments regarding compilation.

Note that NASA's ephemerides are nolonger distributed with the library on
GitHub due to the file size. Additionally, the GSHHG data for plotting maps
and determining if a point is on land is not distributed with the library.
See below for where to download these free datasets.

SOME OF THE CHANGES SINCE VERSION 6.0:
- Bug fixes, notably the CDF in the chi squared distribution now works over
  a much wider range of values of the noncentrality parameter. Previously,
  values as low as 500 were problematic.
- All of the detection statistics classes and related functions have been
  redone to handle a more common normalization and have better comments.
  One can now choose the normalization to use. Also, SwerlingISqLawD now
  has new methods for the first and second derivative of the PDF with
  respect to the average SNR.
- extractSubPlots To let one to pull sub plots out of figures in Matlab
  that have already been plotted with suplots.
- Added functions for displaying maps and for determining whether a point
  is on land or not. For example, run the sample code in the comments to
  pointIsOnLand. Note that these functions rely on the data from the GSHHG,
  which is not provided, but which can be downloaded from
  https://www.ngdc.noaa.gov/mgg/shorelines/ . The binary archive
  gshhg-bin.zip should be placed in Misc/Maps/data.
- Added an example to the comments of expRateLikeGammaConjUpdate on how to
  use the function.
- The inputs and outputs of pseudoInverseDerivative have been reformulated
  to be more useful.
- Added the Erlang distribution in ErlangD and a distribution for a
  weighted sum of noncentral chi squared random variables plus an optional
  normal in ChiSquaredSumD.
- New Lorentz transformation functions along with gradients and Hessians.
- New functions related to simple refraction models. Specifically,
  atmosExpDecayConst4Refrac and approxRefractivity.
- Combinatorial functions for involutions: genAllInvolutions,
  numPureInvolutions, numFixedPointInvolutions, isPureInvolution, and 
  genAllPureInvolutions.
- The C implementation of clipPolygonSH2D was buggy and was removed.
- In all, 55 new functions (not all in unique files) added and one removed.

COMPILED CODE:

The compilation of the library has been tested under Matlab2024a under
Windows 10 using minGW64 and Microsoft Visual C++ 2022. The code will
probably compile under Mac OS X and Linux. Precompiled code is not
distributed with the library.

EXTERNAL SOLVERS:

Almost all of the function in the library work without any external
resources. However, the folder
./Mathematical_Functions/Polynomials/Generic_Multivariate_Polynomials/
contains the function solvePolySysWithExtProg. This is an interface to
external solvers for simultaneous multivariate polynomial systems. The
solvers are external programs and are not included. Though the functions
polyRootsMultiDim and solveQuadBivarEq can solver certain types of
simultaneous multivariate polynomials, they are limited. For more general
systems, an external solver is the best choice. The functions
polyMeasConvert, polyMeasConvertAsync, and DopplerOnlyInit6D are located in
subfolders entitled 'Uses_External_Solver" and use the function
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
https://www.macaulay2.com/

For the external solvers to work, they need to be added to the search path
used in Matlab so that the phs/bertini/m2 commands can be called using the
system command in Matlab.

DATA SOURCES:

If the library was downloaded without associated data files for solar
system ephemerides as well as magnetic, gravitational and terrain models,
then the relevant data files will have to be downloaded The original 
sources of the data are:

1) Earth2014 terrain model:
https://geodesy.curtin.edu.au/research/harmonic-topography/
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

4) The IGRF13 magnetic field coefficients:
http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
The file http://www.ngdc.noaa.gov/IAGA/vmod/igrf13coeffs.txt can be
downloaded and placed into the folder ./Magnetism/data .

5) The WMM2020 magnetic field coefficients:
https://www.ngdc.noaa.gov/geomag/WMM/
The coefficient file WMM2020COF.zip shown be downloaded, unzipped and the
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

10) NASA JPL's DE440t planetary/solar/lunar ephemerides
https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440t/
The file linux_p1550p2650.440t should be downloaded and placed in 
./Astronomical Code/data
for use by the readJPLEphem function. The readJPLEphem function can read
other NASA ephemerides, but the DE440t set are the default.

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

15) A 5400X2700 image of the Earth. The August, Blue Marble Next Generation
    with Topography and Bathymetry from NASA:
https://visibleearth.nasa.gov/images/73776/august-blue-marble-next-generation-w-topography-and-bathymetry
Name the file NASABlueMarble.jpg and place the file in Misc/data . This is
used as the default map to show on the Earth given by the
plotMapOnEllipsoid function.

16) The Global Self-consistent, Hierarchical, High-resolution Geography
    Database
https://www.ngdc.noaa.gov/mgg/shorelines/
Download the zipped binary formatted archive. The file is
gshhg-bin-2.3.7.zip. It should be renamed to gshhg-bin.zip and be placed in
Misc/Maps/data. Note that the claims a weak copyleft license.

September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.<br>
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

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


