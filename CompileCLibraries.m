%%COMPILECLIBRARIES Compile all of the functions implemented in C++ that
%                   are to be called from Matlab. If a C or C++ mex
%                   function has the same name as a file implemented in
%                   Matlab and is compiled, then Matlab will execute the
%                   compiled mex function rather than the (usually slower)
%                   Matlab code.
%
%The file has been sucessfully run on Mac OS X 10.11 and under Windows 10
%with the Microsoft Visual Studio installed. If a compiler used under Windows is
%something other than the one installed with the Microsoft SDK or with
%Microsoft Visual Studio, then the link command to create static libraries
%for the IAU SOFA library and LibAIS might need to be changed. Under Mac OS
%X/ *NIX, make files suffice to compile the libraries. Under Windows, this
%file serves as a makefile. The proper shell environment variables under
%Windows are taken from those used to setup Matlab's mex commmand.
%
%The option '-U__STDC_UTF_16__' is included in all of the mex commands to
%get rid of an error that can occur when using older versions of Matlab
%under Mac OS X.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

function CompileCLibraries()

%Get the path to this file. It should be root level in the tracker library.
ScriptPath=mfilename('fullpath');
ScriptFolder=fileparts(ScriptPath);

%Save the old working directory and switch to the location of this file.
curDir=pwd;
cd(ScriptFolder)

%Compile optimization code
%Compile the SCS library
cd('./3rd_Party_Libraries/scs-1.2.6/matlab')
make_scs

%Move the compiled code into the proper location.
if(isunix()||ismac())%*NIX/ Mac OS X
    ext=mexext();
    system(['mv ./scs_indirect.',ext,' ../../../0_Compiled_Code/scs_indirect.',ext])
    system(['mv ./scs_version.',ext,' ../../../0_Compiled_Code/scs_version.',ext]);
else%Windows
    ext=mexext();
    system(['move scs_indirect.',ext,' ..\..\..\0_Compiled_Code\scs_indirect.',ext]);
    system(['move scs_version.',ext,' ..\..\..\0_Compiled_Code\scs_version.',ext]);
end

cd(ScriptFolder)

%Compile divRectOpt
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/direct','./Mathematical Functions/Continuous Optimization/divRectOpt.c','./3rd_Party_Libraries/direct/direct_wrap.c','./3rd_Party_Libraries/direct/DIRect.c','./3rd_Party_Libraries/direct/DIRserial.c','./3rd_Party_Libraries/direct/DIRsubrout.c');

%Compile lineSearch
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/liblbfgs/include','-I./3rd_Party_Libraries/liblbfgs/lib','./Mathematical Functions/Continuous Optimization/lineSearch.c','./3rd_Party_Libraries/liblbfgs/lib/lbfgs.c');

%Compile quasiNewtonLBFGS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/liblbfgs/include','-I./3rd_Party_Libraries/liblbfgs/lib','./Mathematical Functions/Continuous Optimization/quasiNewtonLBFGS.c','./3rd_Party_Libraries/liblbfgs/lib/lbfgs.c');

%Compile misc code
%Compile code concerning the rounding mode
mex('-v','CFLAGS="$CFLAGS -frounding-math -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Misc/setProcRoundingMode.c');
%Compile index2NDim
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Misc/index2NDim.cpp');
%Compile nDim2Index
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Misc/nDim2Index.cpp');

%Compile navigation code
%Compile indirectGeodeticProb
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/GeographicLib-1.47/include','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Navigation/indirectGeodeticProb.cpp','./Mathematical Functions/Shared C++ Code/wrapRangeCPP.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/Geodesic.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/GeodesicExact.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/GeodesicLine.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/GeodesicLineExact.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/EllipticFunction.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/Math.cpp','./3rd_Party_Libraries/GeographicLib-1.47/src/GeodesicExactC4.cpp');

%Compile general coordinate system code.
%Compile spher2Cart
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/spher2Cart.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp');
%Compile Cart2Sphere
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Cart2Sphere.cpp','./Coordinate Systems/Shared C++ Code/Cart2SphereCPP.cpp');
%Compile ruv2Cart
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/ruv2Cart.cpp','./Coordinate Systems/Shared C++ Code/ruv2CartCPP.cpp');
%Compile Cart2Ruv
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Cart2Ruv.cpp','./Coordinate Systems/Shared C++ Code/Cart2RuvCPP.cpp');
%Compile getRangeRate
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Measurement Components/getRangeRate.cpp','./Coordinate Systems/Shared C++ Code/getRangeRateCPP.cpp');
%Compile state2RuvRR
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/State Conversion/state2RuvRR.cpp','./Coordinate Systems/Shared C++ Code/getRangeRateCPP.cpp','./Coordinate Systems/Shared C++ Code/Cart2RuvCPP.cpp');
%Compile state2SpherRR
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/State Conversion/state2SpherRR.cpp','./Coordinate Systems/Shared C++ Code/getRangeRateCPP.cpp','./Coordinate Systems/Shared C++ Code/Cart2SphereCPP.cpp');
%Compile getENUAxes
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/getENUAxes.cpp','./Coordinate Systems/Shared C++ Code/getENUAxesCPP.cpp');
%Compile getEllipsHarmAxes
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/getEllipsHarmAxes.cpp','./Coordinate Systems/Shared C++ Code/getEllipsHarmAxesCPP.cpp');
%Compile Cart2EllipsHarmon
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Cart2EllipsHarmon.cpp','./Coordinate Systems/Shared C++ Code/Cart2EllipsHarmonCPP.cpp');
%Compile Cart2Ellipse
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Coordinate Systems/Cart2Ellipse.c')

%Compile coordinate system gradient code
%Compile rangeGradient
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Jacobians/Component Gradients/rangeGradient.cpp','./Coordinate Systems/Shared C++ Code/rangeGradientCPP.cpp');
%Compile spherAngGradient
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Jacobians/Component Gradients/spherAngGradient.cpp','./Coordinate Systems/Shared C++ Code/spherAngGradientCPP.cpp');
%Compile calcSpherConvJacob
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Jacobians/Converted Jacobians/calcSpherConvJacob.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeGradientCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngGradientCPP.cpp');
%Compile calcSpherInvJacob
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Jacobians/calcSpherInvJacob.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvJacobCPP.cpp');

%Compile coordinate system Hessian code
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/Component Hessians/rangeHessian.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/Component Hessians/spherAngHessian.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/calcSpherInvHessian.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvHessianCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/calcSpherHessian.cpp','./Coordinate Systems/Shared C++ Code/calcSpherHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/Converted Hessians/calcSpherConvHessian.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp');

%Compile the relativity code
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Relativity/Shared C Code/','./Coordinate Systems/Relativity/relVecAdd.c','./Coordinate Systems/Relativity/Shared C Code/relVecAddC.c');

%Compile the magnetic and gravitational code.
%Compile NALegendreCosRat
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/Polynomials/NALegendreCosRat.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp');
% %Compile LegendreCos
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/Polynomials/LegendreCos.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp');
%Compile normHelmholtz
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/Polynomials/normHelmholtz.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp');
%Compile spherHarmonicEval
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','-I./Coordinate Systems/Shared C++ Code/','./Mathematical Functions/Spherical Harmonics/spherHarmonicEval.cpp','./Mathematical Functions/Shared C++ Code/spherHarmonicEvalCPP.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeGradientCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngGradientCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvHessianCPP.cpp');
%Compile spherHarmonicSetEval
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','-I./Coordinate Systems/Shared C++ Code/','./Mathematical Functions/Spherical Harmonics/spherHarmonicSetEval.cpp','./Mathematical Functions/Shared C++ Code/spherHarmonicSetEvalCPP.cpp','./Mathematical Functions/Shared C++ Code/NALegendreCosRatCPP.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeGradientCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngGradientCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvJacobCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherConvHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvHessianCPP.cpp');
%Compile spherHarmonicCov
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','-I./Coordinate Systems/Shared C++ Code/','./Mathematical Functions/Spherical Harmonics/spherHarmonicCov.cpp','./Mathematical Functions/Shared C++ Code/spherHarmonicCovCPP.cpp','./Mathematical Functions/Shared C++ Code/normHelmholtzCPP.cpp','./Coordinate Systems/Shared C++ Code/spher2CartCPP.cpp');

%Compile the 2D assignment algorithms
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Assignment Algorithms/Association Probabilities and Specific Updates/calc2DAssignmentProbs.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/permCPP.cpp');
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Sample Code/2D Assignment/assign2DByCol.c');
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Assignment Algorithms/2D Assignment/assign2D.c');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/Shared C++ Code/','./Sample Code/2D Assignment/assign2DAlt.cpp','./Assignment Algorithms/Shared C++ Code/ShortestPathCPP.cpp');

%Compile the k-best 2D assignment algorithm
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/Shared C++ Code/','./Assignment Algorithms/k-Best 2D Assignment/kBest2DAssign.cpp','./Assignment Algorithms/Shared C++ Code/ShortestPathCPP.cpp');

%Compile the containers
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Container Classes/metricTreeCPPInt.cpp','./Container Classes/Shared C++ Code/metricTreeCPP.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Container Classes/kdTreeCPPInt.cpp','./Container Classes/Shared C++ Code/kdTreeCPP.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp');

%Compile the mathematical functions
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/Geometry/turnOrientation.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/exactSignOfSum.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/pointIsInPolygon.cpp','./Mathematical Functions/Geometry/Shared C++ Code/pointIsInPolygonCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/twoLineIntersectionPoint2D.cpp','./Mathematical Functions/Geometry/Shared C++ Code/twoLineIntersectionPoint2DCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/signedPolygonArea.cpp','./Mathematical Functions/Geometry/Shared C++ Code/signedPolygonAreaCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C++ Code/','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/clipPolygonSH2D.cpp','./Mathematical Functions/Geometry/Shared C++ Code/twoLineIntersectionPoint2DCPP.cpp','./Mathematical Functions/Geometry/Shared C++ Code/signedPolygonAreaCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/perm.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/permCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/getNextCombo.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/getNextGrayCode.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C++ Code/','-I./Container Classes/Shared C++ Code/','./Mathematical Functions/findFirstMax.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp')
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C Code/','./Mathematical Functions/Searching/binSearchDoubles.c','./Mathematical Functions/Shared C Code/binSearchC.c')
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/Shared C++ Code/','-I./Mathematical Functions/MMOSPAApprox/Shared C++ Code/','./Mathematical Functions/MMOSPAApprox/MMOSPAApprox.cpp','./Mathematical Functions/MMOSPAApprox/Shared C++ Code/MMOSPAApproxCPP.cpp','./Assignment Algorithms/Shared C++ Code/ShortestPathCPP.cpp');
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/wrapRange.cpp','./Mathematical Functions/Shared C++ Code/wrapRangeCPP.cpp')
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Mathematical Functions/Basic Matrix Operations/kronSym.c')

%Functions using LAPACK
lapackInclude=['-I',fullfile(matlabroot,'extern','include')];
if(ispc())
    lapackLib=fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwlapack.lib');
else
    lapackLib='-lmwlapack';
end
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/',lapackInclude,'-I./','./Mathematical Functions/Basic Matrix Operations/tria.c',lapackLib)

%If compiling under Windows, the compile environment must be set up so
%that external libraries can be compiled and linked. The settings that
%Matlab has for the compiler should generally work. If compiling under
%other platforms it is usually not necessary to do anything but add the
%Matlab binaries to the search path so that other scripts can call the mex
%command if they want to compile mex files externally.
if(ispc())
    cc=mex.getCompilerConfigurations();
    compDet=cc.Details;
    %Separate out each command that sets the environment.
    setCommands=strsplit(compDet.SetEnv,'\n');
    numCommands=length(setCommands);
    %While setting the environment, the old values will be saved and the
    %restored befoer this function terminates.
    numSet=0;
    envVars={};
    origEnvVals={};
    for curCommand=1:numCommands
        command=setCommands{curCommand};
        if(~isempty(command))
            %The commands are of the format 'set COMPFLAGS=/c....' However,
            %to set the environment in matlab, one must use the setenv
            %function. This means splitting the command into parts. For
            %example, here, the COMPFLAGS would be the first argument to
            %the setenv command and the string after the first equals sign
            %would be the second argument.
            splitIdx=strfind(command,'=');
            commandString=command((splitIdx+1):end);

            %The first entry contains the 'set' word. Get rid of it and
            %find the parameter being set.
            temp=strsplit(command(1:(splitIdx-1)),' ');
            %The name of the thing being set, like PATH or INCLUDE. We are
            %using end instead of 2, because the strsplit might have slip
            %out whitespace before the word 'set'.
            envVar=temp{end};

            %We want the string of the thing being set. Matlab puts the
            %previous environmental variables at the end of the string.
            %However, that will not work here, because we have to use the
            %getenv command. The previous environmental variable is called
            %with a % sign in front like %PATH%, so we can split on that
            %and just take the first part.
            settingString=strsplit(commandString,'%');
            envVars{numSet+1}=envVar;
            origEnvVals{numSet+1}=getenv(envVar);

            %The mex compile program must be in the path for anyting to
            %compile that wishes to externally call the matlab mex
            %compiler.
            if(strcmp(envVar,'PATH')==1||strcmp(envVar,'Path')==1||strcmp(envVar,'path')==1)
                setenv(envVar,[settingString{1},';',matlabroot,'/bin'])
            else
                setenv(envVar,settingString{1})
            end
            numSet=numSet+1;
        end
    end
else
    numSet=1;
    envVars{1}='PATH';
    origEnvVals{1}=getenv('PATH');
    setenv('PATH',[origEnvVals{1},':',matlabroot,'/bin'])
end

%%%Compile the SOFA library 
%Switch to the SOFA code directory
cd ./3rd_Party_Libraries/sofa/src
%Run commands on the command line to build the library and place it in a
%known location.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('make');
    system('mv ./libsofa_c.a ../../../0_Compiled_Code/libsofa_c.a');
    system('make clean');
    linkCommands{1}='./0_Compiled_Code/libsofa_c.a';
elseif(ispc())
    %Create object files for all of the C++ files.
    cc=[getenv('COMPILER'),' ',getenv('COMPFLAGS'),' ','-U__STDC_UTF_16__ '];
    system(['for %f in (*.c) do ',cc,'%f']) 
    %Store the names of the object files
    system('dir /b *.obj > temp.lst');
    %Get the name of the linker
    ll=[getenv('LINKER'),' '];
    %Combine the object files into a library. it is assumed that the linker
    %takes options akin to the link function in Microsoft's SDK, so that we
    %can tell it to make things into a library rather than trying to
    %compile it into a program.
    system([ll,'-lib /nologo /out:libsofa_c.lib @temp.lst']);
    %Move the library into the folder with the compiled code.
    system('move libsofa_c.lib ..\..\..\0_Compiled_Code\libsofa_c.lib');
    %Get rid of the object files
    system('del *.obj');
    %Get rid of the list of names of the object files.
    system('del temp.lst');
    linkCommands{1}='-llibsofa_c';
    linkCommands{2}='-L./0_Compiled_Code';
else
    error('Cannot Compile SOFA library for unknown operating system')
end

cd(ScriptFolder)

%Compile astronomical functions that use the SOFA code
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Astronomical Code/changeEpoch.c',linkCommands{:})
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./Coordinate Systems/Shared C++ Code/','-I./','./Astronomical Code/starCat2Obs.cpp','./Coordinate Systems/Shared C++ Code/getENUAxesCPP.cpp',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./Coordinate Systems/Relativity/Shared C Code/','-I./','./Astronomical Code/aberrCorr.c','./Coordinate Systems/Relativity/Shared C Code/relVecAddC.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Astronomical Code/lightDeflectCorr.c',linkCommands{:})

%Compile the atmospheric code
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/nrlmsise-00.20151122/','-I./','./Atmospheric Models/NRLMSISE00GasTemp.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00_data.c')
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/nrlmsise-00.20151122/','-I./','./Atmospheric Models/NRLMSISE00Alt4Pres.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00_data.c')
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Atmospheric Models/simpAstroRefParam.c',linkCommands{:})

%%Compile the coordinate transforms that use the SOFA code.
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2ITRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2GCRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2TIRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2GCRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TEME2ITRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2TEME.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TOD2GCRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2TOD.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/MOD2GCRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2MOD.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/J2000F2ICRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2J2000F.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2ITRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2TIRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2CIRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/CIRS2GCRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/CIRS2TIRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2CIRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/G2ICRS.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2G.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2Ecliptic.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ecliptic2ICRS.c',linkCommands{:})

mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/sofa/src/','./Coordinate Systems/ellips2Cart.c',linkCommands{:})

%Compile the time conversion functions that use the SOFA code.
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/UTC2Cal.c','Coordinate Systems/Time/Shared C Code/UTC2CalC.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/UTC2DayCount.c','Coordinate Systems/Time/Shared C Code/UTC2DayCountC.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/Cal2UTC.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/UTC2TAI.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TAI2TT.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TCG2TT.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2TAI.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2TCG.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2GMST.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2GAST.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TAI2UTC.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/BesselEpoch2TDB.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TDB2BesselEpoch.c',linkCommands{:})

mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TDB2TCB.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TCB2TDB.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/cumLeapSec.c',linkCommands{:})

mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/JulDate2JulEpoch.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/JulEpoch2JulDate.c',linkCommands{:})

mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2UT1.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TAI2UT1.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2TCB.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2TDB.c',linkCommands{:})
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TDB2TT.c',linkCommands{:})

%%Compile other astronomical code
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/sofa/src/','./Astronomical Code/approxSolarSysVec.c',linkCommands{:})
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/SGP4/cpp/','-I./','./Astronomical Code/propagateOrbitSGP4.cpp','./3rd_Party_Libraries/SGP4/cpp/sgp4unit.cpp')

%Get rid of the library as nothing else is going to link to it.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('rm ./0_Compiled_Code/libsofa_c.a');
else%Windows
    system('del .\0_Compiled_Code\libsofa_c.lib');
end

%%Compile the AIS library
cd ./3rd_Party_Libraries/libais-0048ceb/src/libais/
%Run commands on the command line to build the library and place it in a
%known location.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('make -f Makefile-custom libais.a');
    system('mv ./libais.a ../../../../0_Compiled_Code/libais.a');
    system('make -f Makefile-custom clean');
    linkCommands{1}='./0_Compiled_Code/libais.a';
elseif(ispc())%Windows    
    %Create object files for all of the C++ files that are used.
    cc=[getenv('COMPILER'),' ',getenv('COMPFLAGS'),' ','-U__STDC_UTF_16__ '];
    system([cc,'ais.cpp']);
    system([cc, 'ais1_2_3.cpp']);
    system([cc, 'ais4_11.cpp']);
    system([cc, 'ais5.cpp']);
    system([cc, 'ais6.cpp']);
    system([cc, 'ais7_13.cpp']);
    system([cc, 'ais8.cpp']);
    system([cc, 'ais8_1_22.cpp']);
    system([cc, 'ais8_1_26.cpp']);
    system([cc, 'ais8_200.cpp']);
    system([cc, 'ais8_366.cpp']);
    system([cc, 'ais8_366_22.cpp']);
    system([cc, 'ais8_367.cpp']);
    system([cc, 'ais9.cpp']);
    system([cc, 'ais10.cpp']);
    system([cc, 'ais12.cpp']);
    system([cc, 'ais14.cpp']);
    system([cc, 'ais15.cpp']);
    system([cc, 'ais16.cpp']);
    system([cc, 'ais17.cpp']);
    system([cc, 'ais18.cpp']);
    system([cc, 'ais19.cpp']);
    system([cc, 'ais20.cpp']);
    system([cc, 'ais21.cpp']);
    system([cc, 'ais22.cpp']);
    system([cc, 'ais23.cpp']);
    system([cc, 'ais24.cpp']);
    system([cc, 'ais25.cpp']);
    system([cc, 'ais26.cpp']);
    system([cc, 'ais27.cpp']);
    system([cc, 'decode_body.cpp']);
    system([cc, 'vdm.cpp']);
    
    %Store the names of the object files
    system('dir /b *.obj > temp.lst');
    %Get the name of the linker
    ll=[getenv('LINKER'),' '];
    %Combine the object files into a library. it is assumed that the linker
    %takes options akin to the link function in Microsoft's SDK, so that we
    %can tell it to make things into a library rather than trying to
    %compile it into a program.
    system([ll,'-lib /nologo /out:libais.lib @temp.lst']);
    %Move the library into the folder with the compiled code.
    system('move libais.lib ..\..\..\..\0_Compiled_Code\libais.lib');
    %Get rid of the object files
    system('del *.obj');
    %Get rid of the list of names of the object files.
    system('del temp.lst');
    linkCommands{1}='-llibais';
    linkCommands{2}='-L./0_Compiled_Code';
else
    error('Cannot Compile AIS library for unknown operating system')
end
cd(ScriptFolder)

%Compile the code to process AIS messages.
mex('-v','CXXFLAGS="$CXXFLAGS -std=c++11"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/libais-0048ceb/src/libais/','-I./','-I./Transponders/Shared C++ Code/','./Transponders/decodeAISString.cpp','./Transponders/Shared C++ Code/AISFuncs.cpp',linkCommands{:})
mex('-v','CXXFLAGS="$CXXFLAGS -std=c++11"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/libais-0048ceb/src/libais/','-I./','-I./Transponders/Shared C++ Code/','./Transponders/decodeAISPosReports2Mat.cpp','./Transponders/Shared C++ Code/AISFuncs.cpp',linkCommands{:})

%Get rid of the library as nothing else is going to link to it.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('rm ./0_Compiled_Code/libais.a');
else%Windows
    system('del .\0_Compiled_Code\libais.lib');
end

%Restore the environment values back to what they originally were.
for curSet=1:numSet
    setenv(envVars{numSet},origEnvVals{numSet});
end

%Restore the old working directory.
cd(curDir);
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
