#!/bin/csh
#
#   Mac OS X Intel cc version.
#
#   This script builds the shared object library for the Matlab-CSPICE
#   interface (Mice). It assumes that it is executed from one of the
#   "product" directories in a tree that looks like the one displayed
#   below:
#
#                      package
#                         |
#                         |
#       +------+------+------+------+------+
#       |      |      |      |      |      |
#     data    doc    etc    exe    lib    src
#                                          |
#                                          |
#                         +----------+----------+------- ... ------+
#                         |          |          |                  |
#                     product_1  product_2  product_3    ...   product_n
#
#   Here's the basic strategy:
#
#     1)  Compile all of the .c files in the current directory to
#         LOCALLIB. The "mex" script include with the Matlab
#         distribution initiates the compile using the options defined
#         in the local mexopts.sh file. The "mex" is assumed to in
#         one of the directories in the current executable path.
#
#     2)  Move LOCALLIB to the 'lib' directory, set the library to the
#         correct name.
#
#   Change History:
#   ===============
#
#   Version 1.0.0,  Feb. 14, 2008, Ed Wright & Boris Semenov
#

#
# Check if "mex" is in the current executable path. If not, complain
# and stop.
#
which mex >& /dev/null

if ( $status != 0 ) then
   echo " "
   echo "      The 'mex' script include with the Matlab is not in the "
   echo "      current path:"
   echo " "
   echo $PATH
   echo " "
   echo "      As a result, I am unable to build this product."
   echo " "
   exit 1;
endif

#
# Set the names of the local and final libraries.
#
set LOCALLIB = "mice.mexmaci64"
set LIBRARY  = "../../lib/$LOCALLIB"

#
# Set compile options.
#
set MEXOPTIONS = "-Dunix"

#
# Compile all .c modules in the current directory into the local
# library.
#
mex $MEXOPTIONS -output $LOCALLIB  -I../../include/ *.c ../../lib/*.a

#
# Move local library to it final location.
#
\mv $LOCALLIB $LIBRARY

#
# All done.
#
exit 0
