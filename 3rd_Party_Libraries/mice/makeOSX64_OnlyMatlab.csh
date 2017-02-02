\echo This script builds the SPICE delivery
\echo for the mice package of the toolkit.
\echo This does not build any of the helper/
\echo sample programs.
\echo
\echo The script must be executed from the
\echo mice directory.
\echo
cd src
\echo
\echo Creating cspice
\echo
cd cspice
chmod u+x mkprodctOSX.csh; ./mkprodctOSX.csh
cd ..
\echo
\echo Creating csupport
\echo
cd csupport
chmod u+x mkprodctOSX.csh; ./mkprodctOSX.csh
cd ..
\echo
\echo Creating mice
\echo
cd mice
chmod u+x mkprodctOSX.csh; ./mkprodctOSX.csh
cd ..
cd ..
\echo Toolkit Build Complete
