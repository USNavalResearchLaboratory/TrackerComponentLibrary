#-Abstract
#
#   The NAIF modified version of the mexopts.sh file provided with
#   MATLAB installations.
#   
#-Disclaimer
#
#   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
#   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
#   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
#   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
#   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
#   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
#   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
#   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
#   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
#   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
#
#   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
#   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
#   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
#   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
#   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
#   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
#   KNOW OF THE POSSIBILITY.
#
#   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
#   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
#   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
#   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
#
#-I/O
#
#   None
#
#-Examples
#
#   None.
#
#-Particulars
#
#   Mathworks provides a mexopts.sh file to define the platform
#   specific options for compiling/linking MEX libraries. NAIF provides
#   this version of the file with edits to build the Mice interface for
#   platforms:
#
#      OS X Intel GCC 32-bit -> maci
#      OS X PPC GCC          -> mac
#      Linux GCC 32-bit      -> glnx86
#      Linux GCC 64-bit      -> glnxa64
#      Solaris cc 32-bit     -> sol2
#      Solaris cc 64-bit     -> sol64
#      OS X Intel GCC 64-bit -> maci64
#
#   -fexceptions
#
#           Enable exception handling.  Generates extra code needed to propa-
#           gate exceptions.  For some targets, this implies GCC will generate
#           frame unwind information for all functions, which can produce sig-
#           nificant data size overhead, although it does not affect execution.
#           If you do not specify this option, GCC will enable it by default
#           for languages like C++ which normally require exception handling,
#           and disable it for languages like C that do not normally require
#           it.  However, you may need to enable this option when compiling C
#           code that needs to interoperate properly with exception handlers
#           written in C++.  You may also wish to disable this option if you
#           are compiling older C++ programs that don't use exception handling.
#
#
#   mexopts.sh
#
#          Shell script for configuring MEX-file creation script,
#          mex.  These options were tested with the specified compiler.
#
#   usage: Do not call this file directly; it is sourced by the
#          mex shell script.  Modify only if you don't like the
#          defaults after running mex.  No spaces are allowed
#          around the '=' in the variable assignment.
#
#-Required Reading
#
#   None.
#
#-Version
#
#   -Mice 1.2.0 16-FEB-2010 (EDW)
#
#      Added sol64 section for Solaris cc 64-bit, copied from default
#      Matlab mexopts.sh.
#
#      Added sol2 section for Solaris cc 32bit, copied from default
#      Matlab mexopts.sh.
#
#      Added maci64 section for OS X cc 64-bit, copied from default
#      Matlab mexopts.sh.
#
#   -Mice 1.1.0 04-FEB-2009 (EDW)
#
#      Added glnxa64 section, copied from default Matlab mexopts.sh.
#
#   -Mice 1.0.1 10-JAN-2008 (EDW)
#
#      Added -O2 to Linux compile paramters.
#
#   -Mice 1.0.0 10-JAN-2008 (EDW)
#
#&


#
# Define the relative location of the CSPICE root directory and the 
# additional compiler warnings
#
SPICE='../..'
GCC_WARNING=' -ansi -Wall -Wundef -Wpointer-arith -Wcast-align -Wsign-compare ' 



#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files via the system ANSI compiler
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision: 1.85 $  $Date: 2002/06/10 18:56:23 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmatlb -lmat -lmwservices -lut -lm"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lm"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,--rpath-link,$TMW_ROOT/extern/lib/$Arch,--rpath-link,$TMW_ROOT/bin/$Arch"
#           gcc -v
#           gcc version 2.95.2 19991024 (release)

            CC='gcc'
            CFLAGS="-fPIC -O2 -ansi -D_GNU_SOURCE -pthread -fexceptions  -I$SPICE/include $GCC_WARNING"
            CLIBS=" $SPICE/lib/cspice.a $RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#           
#           g++ -v
#           gcc version 2.95.2 19991024 (release)
            CXX='g++'
#           Add -fhandle-exceptions to CXXFLAGS to support exception handling
            CXXFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread'
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#           g77 -v -xf77-version 
#           g77 version 2.95.2 19991024 (release) 
#           (from FSF-g77 version 0.5.25 19991024 (release))
#           NOTE: g77 is not thread safe
            FC='g77'
            FFLAGS='-fPIC'
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE -fexceptions'
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CFLAGS="$CFLAGS -O2  -I$SPICE/include $GCC_WARNING"
            CLIBS="$RPATH $MLIBS $SPICE/lib/cspice.a  -lm -lstdc++"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexa64'
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS="-fno-common -no-cpp-precomp -O2 -fexceptions -I$SPICE/include $GCC_WARNING"

            CLIBS=" $SPICE/lib/cspice.a $MLIBS "

            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'

            if [ -f /usr/bin/g++2 ]; then 
                CXX=g++2
            else
                CXX=c++
            fi
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='f77'
            FFLAGS='-f -N15 -N11 -s -Q51 -W'
            ABSOFTLIBDIR=`which $FC | sed -n -e '1s|bin/'$FC'|lib|p'`
            FLIBS="-L$ABSOFTLIBDIR -lfio -lf77math"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS='nmedit -s $TMW_ROOT/extern/lib/$Arch/$MAPFILE $mex_file'
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS="-fno-common -no-cpp-precomp -O2 -arch i386 -fexceptions -I$SPICE/include $GCC_WARNING"

            CLIBS=" $SPICE/lib/cspice.a  $MLIBS "

            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'

            if [ -f /usr/bin/g++2 ]; then 
                CXX=g++2
            else
                CXX=c++
            fi
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch i386'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='g95'
            FFLAGS=''
            FLIBS="$MLIBS"
            FOPTIMFLAGS='-O '
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci'
            LDFLAGS="-bundle -arch i386 -Wl,-flat_namespace -undefined suppress  -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            # CkeyName: Clang
            # CkeyManufacturer: Apple
            # CkeyLanguage: C
            # CkeyVersion:
            # CkeyLinkerName:
            # CkeyLinkerVersion:
            CC='xcrun  -sdk macosx10.9 clang'
## workaround clang defect temporarily use line below           SDKROOT='/Developer/SDKs/MacOSX10.9.sdk'
# compute SDK root on the fly
# target 10.8
            MW_SDKROOT_TMP="find `xcode-select -print-path` -name MacOSX10.9.sdk"
			MW_SDKROOT=`$MW_SDKROOT_TMP`
            MACOSX_DEPLOYMENT_TARGET='10.9'
            ARCHS='x86_64'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -I$SPICE/include -I/usr/include $GCC_WARNING"
            CFLAGS="$CFLAGS  -fexceptions"
            CLIBS=" $SPICE/lib/cspice.a  $MLIBS "
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            # C++keyName: Clang++
            # C++keyManufacturer: Apple
            # C++keyLanguage: C++
            # C++keyVersion:
            # C++keyLinkerName:
            # C++keyLinkerVersion:
            CXX='xcrun  -sdk macosx10.9  clang++'
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -I/usr/include"
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: GNU Fortran
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='gfortran'
            FFLAGS='-fexceptions -m64 -fbackslash'
            FC_LIBDIR=`$FC -print-file-name=libgfortran.dylib 2>&1 | sed -n '1s/\/*libgfortran\.dylib//p'`
            FC_LIBDIR2=`$FC -print-file-name=libgfortranbegin.a 2>&1 | sed -n '1s/\/*libgfortranbegin\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lgfortran -L$FC_LIBDIR2 -lgfortranbegin"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci64'
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
#           cc -V
#           WorkShop Compilers 5.0 98/12/15 C 5.0
            CC='cc'
            CFLAGS='-KPIC -dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt -I$SPICE/include $WARNING'
            CLIBS="$SPICE/lib/cspice.a $MLIBS -lm -lc"
            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-g'
#           
#           CC -V
#           WorkShop Compilers 5.0 98/12/15 C++ 5.0
            CXX='CC -compat=5'
            CCV=`CC -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CXXFLAGS='-KPIC -dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt'
            CXXLIBS="$MLIBS -lm -lCstd -lCrun"
            CXXOPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#           f77 -V
#           WorkShop Compilers 5.0 99/09/16 FORTRAN 77 5.0 patch 107596-03
            FC='f77'
            FFLAGS='-KPIC -dalign -mt'
            FLIBS="$MLIBS -lF77 -lM77 -lsunmath -lm -lcx -lc"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-G -mt -M$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
            CC='cc -m64'
            CFLAGS='-KPIC -dalign -v -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt -I$SPICE/include'
            CLIBS="$SPICE/lib/cspice.a $MLIBS -lm -lc"
            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-xs -g'
#           
            CXX='CC -m64 -compat=5'
            CCV=`CC -m64 -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CXXFLAGS='-KPIC -dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt -library=stlport4,Crun'
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CXXDEBUGFLAGS='-xs -g'
#
            FC='f90 -m64'
            FFLAGS='-KPIC -dalign -mt'
            FLIBS="$MLIBS -lfui -lfsu -lsunmath -lm -lc"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-xs -g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexs64'
            LDFLAGS="-G -mt -M$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-xs -g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac

