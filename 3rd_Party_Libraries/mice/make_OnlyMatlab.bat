rem This script builds the SPICE delivery
rem for the mice package of the toolkit.
rem This does not build any of the helper/
rem sample programs.
rem
rem The script must be executed from the
rem mice directory.
rem
cd src
rem
rem Creating cspice
rem
cd cspice
call mkprodct.bat
cd ..
rem
rem Creating csupport
rem
cd csupport
call mkprodct.bat
cd ..
rem
rem Creating mice
rem
cd mice
call mkprodct.bat
cd ..
cd ..
rem Toolkit Build Complete
