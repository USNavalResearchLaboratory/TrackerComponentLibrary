**Tracker Component Library Release 6.1.1, June 2025**
https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary

A paper describing a number of features of the library is<br>
D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid
Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, no.
5, pp. 18-27, May. 2017.

Those looking for magnetic field synthesis code might want to look at<br>
./Sample Code/Magnetic Models/ .

These are the release notes for the version 6.11 of the Tracker Component
Library. The Tracker Component Library is a collection of Matlab routines
for simulating and tracking targets in various scenarios. Due to the
complexity of the target tracking problem, a great many routines can find
use in other areas including combinatorics, astronomy, and statistics.

Making this code available and open source is in the spirit of the Federal
Source Code Policy, which is OMB M-16-21, which is linked at<br>
https://open.gsa.gov/oss-policy/<br>
https://code.gov/

Those looking to get a quick idea of a very simple end-to-end tracking
algorithm (track initiation, data association, maintenance, and
termination) might want to look at<br>
./Sample Code/Basic Tracking Examples/demo2DIntegratedDataAssociation.m<br>
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

SOME OF THE CHANGES SINCE VERSION 6.1:

- Bug fixes. This includes some in the compiled implementations as well as
  other things. For example, the CDF of the Poisson distribution now
  functions correctly over a much wider range of parameters.
- Functions for solving the Josephus problem and generating Josephus
  permutations.

COMPILED CODE:

The compilation of the library has been tested under Matlab2024a under
Windows 11 using minGW64 and Microsoft Visual C++ 2022. The code will
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

The solvers can be downloaded from:<br>
PHCpack for Mac OS X, Windows and Linux:<br>
http://homepages.math.uic.edu/~jan/download.html<br>
Bertini (version 1.5) for Mac OS X and Linux and Windows, only with Cygwin:<br>
https://bertini.nd.edu/download.html<br>
Macaulay2 for Mac OS X, Linux, and Windows:<br>
https://www.macaulay2.com/

For the external solvers to work, they need to be added to the search path
used in Matlab so that the phs/bertini/m2 commands can be called using the
system command in Matlab.

DATA SOURCES:

If the library was downloaded without associated data files for solar
system ephemerides as well as magnetic, gravitational and terrain models,
then the relevant data files will have to be downloaded The original 
sources of the data are:

1) Earth2014 terrain model:<br>
https://geodesy.curtin.edu.au/research/harmonic-topography/<br>
The files Earth2014.BED2014.degree2160.bshc,
Earth2014.ICE2014.degree2160.bshc, Earth2014.RET2014.degree2160.bshc,
Earth2014.SUR2014.degree2160.bshc, Earth2014.TBI2014.degree2160.bshc should
be placed into ./Terrain/data . Though data files to degree and order
10,800 are also available, they are not currently supported as the function
spherHarmonicEval does not use extended precision arithmetic, which is
required to avoid overflows when using the model.

2) EGM2008 terrain model (The DTM2006.0 model):<br>
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html<br>
The file Coeff_Height_and_Depth_to2190_DTM2006.0 should be placed into 
./Terrain/data . The first time using the getEGM2008TerrainCoeffs function
will be significantly slower than subsequent calls, assuming the full model
is loaded the first time.

3) The EMM2017 magnetic field coefficients:<br>
https://www.ngdc.noaa.gov/geomag/EMM/<br>
The coefficients are included with the software that is provided. The files
having names like EMM2015.COF and EMM2000.COF should all be zipped into an
archive with the name EMM_Coefficients.zip and placed in ./Magnetism/data .
The files with names like EMM2015SV.COF andEMM2000SV.COF should all be
zipped into an archive with the name EMM_Secular_Variations.zip and placed
into ./Magnetism/data . If all of the data for a single year is loaded,
then a .mat file containing the data is placed in  the data folder to speed
up subsequent calls to the function.

4) The IGRF14 magnetic field coefficients:<br>
https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field<br>
The file igrf14coeffs.txt can be downloaded and placed into the folder
./Magnetism/data .

5) The WMM and the WMMHR 2025 magnetic field coefficients for 2025:<br>
https://www.ngdc.noaa.gov/geomag/WMM/<br>
https://www.ncei.noaa.gov/products/world-magnetic-model-high-resolution<br>
The coefficient files WMM2025COF.zip and WMMHR2025COF.zip should be
downloaded, unzipped and the files WMM.COF and WMMHR.COF placed into
./Magnetism/data .

6) EGM2008 gravitational coefficients:<br>
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html<br>
Download the file  EGM2008_to2190_TideFree.gz, which is under "EGM2008 Tide
Free Spherical Harmonic Coefficients". Decompress the file and place it in
./Gravity/data .

7) EGM96 gravitational coefficients:<br>
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html<br>
The file egm96.z should be downloaded and decompressed, resulting in a file
names egm96 without an extension, which should be placed in
./Gravity/data .

8) EGM2008 parameters needed to convert height anomalies to geoid
   undulations:<br>
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html<br>
The "Correction Model" file Zeta-to-N_to2160_egm2008.gz should be
downloaded, ungzipped and placed in ./Gravity/data .

9) EGM96 parameters needed to convert height anomalies to geoid
   undulations:<br>
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html<br>
The file corrcoef.z should be decompressed and placed in ./Gravity/data .

10) NASA JPL's DE440t planetary/solar/lunar ephemerides:<br>
https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440t/<br>
The file linux_p1550p2650.440t should be downloaded and placed in 
./Astronomical Code/data
for use by the readJPLEphem function. The readJPLEphem function can read
other NASA ephemerides, but the DE440t set are the default.

11) NASA JPL's GL0900C gravitational model of the Moon:<br>
http://pds-geosciences.wustl.edu/missions/grail/default.htm
Download the file jggrx_0900c_sha.tab and place it in ./Gravity/data .

12) The FES2004 tide model of the effects of ocean Earth tides:<br>
http://maia.usno.navy.mil/conv2010/convupdt/convupdt_c6.html<br>
The file fes2004_Cnm-Snm.dat should be downloaded and placed in
./Gravity/Tides/data .

14) The Hipparcos 2 Star Catalog:<br>
http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/311<br>
Place the file hip2.dat in ./Astronomical Code/data .

15) A 5400X2700 image of the Earth. The August, Blue Marble Next Generation
    with Topography and Bathymetry from NASA:<br>
https://visibleearth.nasa.gov/images/73776/august-blue-marble-next-generation-w-topography-and-bathymetry<br>
Name the file NASABlueMarble.jpg and place the file in Misc/data . This is
used as the default map to show on the Earth given by the
plotMapOnEllipsoid function.

16) The Global Self-consistent, Hierarchical, High-resolution Geography
    Database:<br>
https://www.soest.hawaii.edu/pwessel/gshhg/<br>
Download the zipped binary formatted archive. The file is
gshhg-bin-2.3.7.zip. It should be renamed to gshhg-bin.zip and be placed in
Misc/Maps/data. Note that the claims a weak copyleft license.

May 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.<br>
(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
