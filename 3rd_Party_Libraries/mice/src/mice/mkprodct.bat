rem
rem  mkprodct.bat for 64-bit Mice, Windows.
rem
rem Version 
rem
rem    1.0.0, May 11, 2010, Ed Wright
rem

rem
rem  Set library name.
rem 

set NAME=mice

rem
rem  Set file names based on the library name.
rem



rem

rem



rem
rem 
rem 

date /t 
time /t 


rem
rem Build the MEX. The mexopts.bat file contains the needed compile and link settings
rem from the build and link against CSPICE. Maintain mexopts.bat file in the source
rem directory with the mkprodct.bat script.
rem
rem The script writes the MEX directly to the "lib" directory.
rem

mex -output ..\..\lib\mice LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib" mice.c  zzmice.c  zzmice_CreateIntScalar.c -I..\..\include ..\..\lib\cspice.lib

