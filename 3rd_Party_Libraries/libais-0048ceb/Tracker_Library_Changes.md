The Tracker Component Library does not use all of the files in aislib.
Unused code was placed into ./Unused Code.zip, ./src/Unused Code.zip, and
./src/libais/Unused Code.zip The executable in the bin folder was removed.

In the file ./src/libais/Makefile-custom, the -fPIC C++ flag was added to
eliminate compile problems on some systems. The -DNDEBUG flasg was also
added to eliminate parts related to debugging. A number of modifications
were made to get rid of warnings when compiling in Visual Studio and MinGW.

The code was originally taken from
https://github.com/schwehr/libais

May 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.