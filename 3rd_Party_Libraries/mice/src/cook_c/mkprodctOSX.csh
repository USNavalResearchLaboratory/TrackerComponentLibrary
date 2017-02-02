#! /bin/csh
#
#   Mac OS X cc version. cook_c version.
#
#   This script is a more or less generic library/executable
#   builder for CSPICE products.  It assumes that it is executed
#   from one of the "product" directories in a tree that looks like
#   the one displayed below:
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
#     1)  Compile all of the .c files in the current directory
#
#     2)  If there are no .pgm files in the current directory this
#         is assumed to be a library source directory.  The name
#         of the library is the same as the name of the product.
#         The library is placed in the "lib" directory in the tree
#         above.  The script is then done.
#
#         If there are .pgm files and there were some .c
#         files compiled the objects are gathered together in the
#         current directory into a library called locallib.a.
#
#     3)  If any *.pgm files exist in the current directory, compile
#         them and add their objects to locallib.a.  Create a C main
#         program file from the uniform CSPICE main program main.x.
#         Compile this main program and link its object with locallib.a,
#         ../../cspice.a and ../../csupport.a. The output
#         executables have an empty extension.  The executables are
#         placed in the "exe" directory in the tree above.
#
#   The environment variable TKCOMPILEOPTIONS containing compile options
#   is optionally set. If it is set prior to executing this script,
#   those options are used. It it is not set, it is set within this
#   script as a local variable.
#
#   References:
#   ===========
#
#   "Unix Power Tools", page 11.02
#      Use the "\" character to unalias a command temporarily.
#
#   "A Practical Guide to the Unix System"
#
#   "The Unix C Shell Field Guide"
#
#   Change History:
#   ===============
#
#   Version 1.1.1  Jul. 30, 2002  Boris Semenov 
#
#      Misc. cleanup: moved "rm locallib.a" to the end of the script, 
#      eliminated duplicate "echo"s, etc. etc.
#
#   Version 1.1.0  Jan. 18, 2002  Ed Wright
#
#      Updated GCC version for Mac OS X CC.
#
#   Version 1.0.0  Feb. 26, 1999  Nat Bachman
#
#

#
#  Choose your compiler.
#
if ( $?TKCOMPILER ) then

    echo " "
    echo "      Using compiler: "
    echo "      $TKCOMPILER"

else

    set TKCOMPILER  =  "cc"
    echo " "
    echo "      Setting default compiler:"
    echo $TKCOMPILER

endif


#
#  What compile options do we want to use? If they were
#  set somewhere else, use those values.  The same goes
#  for link options.
#
if ( $?TKCOMPILEOPTIONS ) then
    echo " "
    echo "      Using compile options: "
    echo "      $TKCOMPILEOPTIONS"
else
#
#  Options:
#
#     -O2                optimize
#
#     -DNON_UNIX_STDIO   Don't assume standard Unix stdio.h
#                        implementation
#
#
    set TKCOMPILEOPTIONS = "-m64 -c -ansi -O2 -DNON_UNIX_STDIO"
    echo " "
    echo "      Setting default compile options:"
    echo "      $TKCOMPILEOPTIONS"
endif

if ( $?TKLINKOPTIONS ) then
    echo " "
    echo "      Using link options: "
    echo "      $TKLINKOPTIONS"
else
    set TKLINKOPTIONS = "-m64 -lm"
    echo " "
    echo "      Setting default link options:"
    echo "      $TKLINKOPTIONS"
endif

echo " "

#
#   Determine a provisional LIBRARY name.
#
foreach item ( `pwd` )
    set LIBRARY = "../../lib/"$item:t
end

#
#  Are there any *.c files that need to be compiled?
#
\ls *.c >& /dev/null

if ( $status == 0 ) then

    foreach SRCFILE ( *.c )
       echo "      Compiling: "   $SRCFILE
       $TKCOMPILER $TKCOMPILEOPTIONS $SRCFILE
    end

endif

echo " "

#
#  If object files exist, we need to create an object library.
#

\ls *.pgm >& /dev/null

if ( $status == 0 ) then
    set LIBRARY = "locallib"
endif

\ls *.o >& /dev/null

if ( $status == 0 ) then

    echo "      Inserting objects in the library $LIBRARY ..."
    ar  crv $LIBRARY.a *.o
    ranlib  $LIBRARY.a
    \rm                *.o
    echo " "

endif

#
#  If there are any main programs in the directory, compile
#  them. If they have their own locallib.a link with it in addition
#  to the default libraries.
#

\ls *.pgm >& /dev/null

if ( $status == 0 ) then

    echo " "

    foreach MAIN ( *.pgm )

       set STEM    = $MAIN:r
       set TARGET  = $STEM.c
       set MAINOBJ = $STEM.o
       set EXECUT = ../../exe/$STEM

       \cp $MAIN $TARGET

       echo "      Compiling and linking: " $MAIN

       if ( -e locallib.a ) then

          $TKCOMPILER    $TKCOMPILEOPTIONS $TARGET
          $TKCOMPILER -o $EXECUT           $MAINOBJ             \
                                           locallib.a           \
                                           ../../lib/csupport.a \
                                           ../../lib/cspice.a   \
                                           $TKLINKOPTIONS

          \rm $TARGET
          \rm $MAINOBJ

       else

          $TKCOMPILER    $TKCOMPILEOPTIONS $TARGET
          $TKCOMPILER -o $EXECUT           $MAINOBJ             \
                                           ../../lib/csupport.a \
                                           ../../lib/cspice.a   \
                                          $TKLINKOPTIONS

          \rm $TARGET
          \rm $MAINOBJ

       endif

    end

endif

#
#  Cleanup.
#

echo " "

\ls *.o >& /dev/null

if ( $status == 0 ) then
    \rm *.o
endif

\ls locallib.a* >& /dev/null

if ( $status == 0 ) then
   \rm locallib.a*
endif


exit 0
