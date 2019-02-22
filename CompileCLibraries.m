%%COMPILECLIBRARIES Compile all of the functions implemented in C or C++
%                   that are to be called from Matlab in the unlimited
%                   distribution portion of the Tracker Component Library.
%                   If a C or C++ mex function has the same name as a file
%                   implemented in Matlab and is compiled, then Matlab will
%                   execute the compiled mex function rather than the
%                   (usually slower) Matlab code if it finds the compiled
%                   version first in its path. Thus, the compiled code is
%                   kept in a 0_Compiled_Code folder at the root level of
%                   the library to ensure that the compiled versions are
%                   found first.
%
%The file has been sucessfully run on Mac OS X 10.12 using XCode
%(clang/llvm) as the compiled and under CentOS Linux using the gcc
%compiler.
%
%The option '-U__STDC_UTF_16__' is included in all of the mex commands to
%get rid of an error that can occur when using older versions of Matlab
%under Mac OS X.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
cd('./3rd_Party_Libraries/scs-2.0.2/scs-matlab-master')
make_scs

%Move the compiled code into the proper location.
if(isunix()||ismac())%*NIX/ Mac OS X
    ext=mexext();
    system(['mv ./scs_indirect.',ext,' ../../../0_Compiled_Code/scs_indirect.',ext]);
    system(['mv ./scs_version.',ext,' ../../../0_Compiled_Code/scs_version.',ext]);
else%Windows
    ext=mexext();
    system(['move scs_indirect.',ext,' ..\..\..\0_Compiled_Code\scs_indirect.',ext]);
    system(['move scs_version.',ext,' ..\..\..\0_Compiled_Code\scs_version.',ext]);
end

cd(ScriptFolder)

%Compile divRectOpt
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/direct','./Mathematical Functions/Continuous Optimization/divRectOpt.c','./3rd_Party_Libraries/direct/direct_wrap.c','./3rd_Party_Libraries/direct/DIRect.c','./3rd_Party_Libraries/direct/DIRserial.c','./3rd_Party_Libraries/direct/DIRsubrout.c');

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
%Compile rangeHessian
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/Component Hessians/rangeHessian.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp');
%Compile spherAngHessian
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/Component Hessians/spherAngHessian.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp');
%Compile calcSpherInvHessian
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/calcSpherInvHessian.cpp','./Coordinate Systems/Shared C++ Code/calcSpherInvHessianCPP.cpp');
%Compile calcSpherHessian
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Coordinate Systems/Shared C++ Code/','./Coordinate Systems/Hessians/calcSpherHessian.cpp','./Coordinate Systems/Shared C++ Code/calcSpherHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/rangeHessianCPP.cpp','./Coordinate Systems/Shared C++ Code/spherAngHessianCPP.cpp');
%Compile calcSpherConvHessian
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
%Compile calc2DAssignmentProbs
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Assignment Algorithms/Association Probabilities and Specific Updates/calc2DAssignmentProbs.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/permCPP.cpp');
%Compile assign2D
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/2D Assignment/Shared C Code/','./Assignment Algorithms/2D Assignment/assign2D.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DC.c');
%Compile assign2DMissedDetect
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/2D Assignment/Shared C Code/','./Assignment Algorithms/2D Assignment/assign2DMissedDetect.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DMissedDetectC.c');
%Compile assign2DFull
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/2D Assignment/Shared C Code/','./Assignment Algorithms/2D Assignment/assign2DFull.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DFullC.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DMissedDetectC.c');
%Compile assign2DByCol
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Sample Code/2D Assignment/assign2DByCol.c');
%Compile assign2DAlt
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/','./Sample Code/2D Assignment/assign2DAlt.cpp','./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/ShortestPathCPP.cpp');
%Compile kBest2DAssign
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/','./Assignment Algorithms/k-Best 2D Assignment/kBest2DAssign.cpp','./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/ShortestPathCPP.cpp');

%Compile the containers
%Compile metricTreeCPPInt
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Container Classes/metricTreeCPPInt.cpp','./Container Classes/Shared C++ Code/metricTreeCPP.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp');
%Compile kdTreeCPPInt
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Container Classes/kdTreeCPPInt.cpp','./Container Classes/Shared C++ Code/kdTreeCPP.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp');

%Compile the mathematical functions
%Compile turnOrientation
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Accurate Arithmetic/Shared C++ Code/','./Mathematical Functions/Geometry/turnOrientation.cpp');
%Compile exactSignOfSum
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Accurate Arithmetic/Shared C++ Code/','./Mathematical Functions/Accurate Arithmetic/exactSignOfSum.cpp');
%Compile pointIsInPolygon
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/pointIsInPolygon.cpp','./Mathematical Functions/Geometry/Shared C++ Code/pointIsInPolygonCPP.cpp');
%Compile twoLineIntersectionPoint2D
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/twoLineIntersectionPoint2D.cpp','./Mathematical Functions/Geometry/Shared C++ Code/twoLineIntersectionPoint2DCPP.cpp');
%Compile signedPolygonArea
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/signedPolygonArea.cpp','./Mathematical Functions/Geometry/Shared C++ Code/signedPolygonAreaCPP.cpp');
%Compile clipPolygonSH2D
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Accurate Arithmetic/Shared C++ Code/','-I./Mathematical Functions/Geometry/Shared C++ Code/','./Mathematical Functions/Geometry/clipPolygonSH2D.cpp','./Mathematical Functions/Geometry/Shared C++ Code/twoLineIntersectionPoint2DCPP.cpp','./Mathematical Functions/Geometry/Shared C++ Code/signedPolygonAreaCPP.cpp');
%Compile perm
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/perm.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/permCPP.cpp');
%Compile getNextCombo
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/getNextCombo.cpp','./Mathematical Functions/Combinatorics/Shared C++ Code/getNextComboCPP.cpp');
%Compile getNextGrayCode
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Combinatorics/Shared C++ Code/','./Mathematical Functions/Combinatorics/getNextGrayCode.cpp');
%Compile findFirstMax
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C++ Code/','-I./Container Classes/Shared C++ Code/','./Mathematical Functions/findFirstMax.cpp','./Mathematical Functions/Shared C++ Code/findFirstMaxCPP.cpp')
%Compile binSearchDoubles
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C Code/','./Mathematical Functions/Searching/binSearchDoubles.c','./Mathematical Functions/Shared C Code/binSearchC.c')
%Compile MMOSPAApprox
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/','-I./Mathematical Functions/MMOSPAApprox/Shared C++ Code/','./Mathematical Functions/MMOSPAApprox/MMOSPAApprox.cpp','./Mathematical Functions/MMOSPAApprox/Shared C++ Code/MMOSPAApproxCPP.cpp','./Assignment Algorithms/k-Best 2D Assignment/Shared C++ Code/ShortestPathCPP.cpp');
%Compile wrapRange
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Container Classes/Shared C++ Code/','-I./Mathematical Functions/Shared C++ Code/','./Mathematical Functions/wrapRange.cpp','./Mathematical Functions/Shared C++ Code/wrapRangeCPP.cpp')
%Compile kronSym
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','./Mathematical Functions/Basic Matrix Operations/kronSym.c')
%Compile heapSortVec
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Shared C Code/','./Mathematical Functions/Sorting/heapSortVec.c','./Mathematical Functions/Shared C Code/heapSortVecC.C')
%Compile permuteMatrix
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Basic Matrix Operations/Shared C Code/','./Mathematical Functions/Basic Matrix Operations/permuteMatrix.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/permuteMatrixC.c')
%Compile minMatOverDim
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Mathematical Functions/Basic Matrix Operations/Shared C Code/','./Mathematical Functions/Basic Matrix Operations/minMatOverDim.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/minMatOverDimC.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/basicMatOps.c');

%Functions using LAPACK
lapackInclude=['-I',fullfile(matlabroot,'extern','include')];
if(ispc())
    lapackLib=fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwlapack.lib');
else
    lapackLib='-lmwlapack';
end
%Compile tria
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/',lapackInclude,'-I./','./Mathematical Functions/Basic Matrix Operations/tria.c',lapackLib)

%Functions using BLAS
blasInclude=['-I',fullfile(matlabroot,'extern','include')];
blasLib='-lmwblas';

%Compile the 3D assignment algorithms.
%Compile assign3D
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/',blasInclude,'-I./','-I./Assignment Algorithms/2D Assignment/Shared C Code/','-I./Assignment Algorithms/3D Assignment/Shared C Code/','-I./Mathematical Functions/Basic Matrix Operations/Shared C Code/','./Assignment Algorithms/3D Assignment/assign3D.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DC.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/basicMatOps.c','./Assignment Algorithms/3D Assignment/Shared C Code/assign3DC.c',blasLib);

%Compile assign3DLB
mex('-v','CFLAGS="$CFLAGS -std=c99 -Wall"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./Assignment Algorithms/2D Assignment/Shared C Code/','-I./Assignment Algorithms/3D Assignment/Shared C Code/','-I./Mathematical Functions/Basic Matrix Operations/Shared C Code/','./Assignment Algorithms/3D Assignment/assign3DLB.c','./Assignment Algorithms/2D Assignment/Shared C Code/assign2DC.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/basicMatOps.c','./Mathematical Functions/Basic Matrix Operations/Shared C Code/minMatOverDimC.c','./Assignment Algorithms/3D Assignment/Shared C Code/assign3DLBC.c');

%%%Compile the SOFA library
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

%Switch to the SOFA code directory
cd ./3rd_Party_Libraries/sofa/src
%Run commands on the command line to build the library and place it in a
%known location.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('make -e CFLAGF=''-c -pedantic -Wall -W -O -fPIC'' CFLAGX=''-pedantic -Wall -W -O -fPIC''');
    system('mv ./libsofa_c.a ../../../0_Compiled_Code/libsofa_c.a');
    system('make clean');
    linkCommands{1}='./0_Compiled_Code/libsofa_c.a';
elseif(ispc())
    compPath=getenv('COMPILER');
    [~,compName]=fileparts(compPath);
    cc=[compPath,' ',getenv('COMPFLAGS'),' ','-U__STDC_UTF_16__ '];
    ll=getenv('LINKER');
    
    if(strcmp(compName,'gcc'))%If here, then minGW is probably being used.
        system(['for %f in (*.c) do ',cc,'%f -c'])
        
        %Store the names of the object files
        system('dir /b *.o > temp.lst');
        
        %Link the object files into an archive.
        system('ar rv ./libsofa_c.a @temp.lst');
        system('ranlib ./libsofa_c.a');
        
        %Move the library into the folder with the compiled code.
        system('move libsofa_c.a ..\..\..\0_Compiled_Code\libsofa_c.a');
        
        %Get rid of the object files
        system('del *.o');
        
        %Get rid of the list of names of the object files.
        system('del temp.lst');
        linkCommands{1}='./0_Compiled_Code/libsofa_c.a';
    else%If here, Microsoft Visual Studio is probably being used.
        %Create object files for all of the C files.
        system(['for %f in (*.c) do ',cc,'%f'])
        
        %Store the names of the object files
        system('dir /b *.obj > temp.lst');
        
        %Combine the object files into a library. it is assumed that the
        %linker takes options akin to the link function in Microsoft Visual
        %Studio, so that we can tell it to make things into a library
        %rather than trying to compile it into a program.
        system([ll,' -lib /nologo /out:libsofa_c.lib @temp.lst']);
    
        %Move the library into the folder with the compiled code.
        system('move libsofa_c.lib ..\..\..\0_Compiled_Code\libsofa_c.lib');
        %Get rid of the object files
        system('del *.obj');
        %Get rid of the list of names of the object files.
        system('del temp.lst');
        linkCommands{1}='-llibsofa_c';
        linkCommands{2}='-L./0_Compiled_Code';
    end
else
    error('Cannot Compile SOFA library for unknown operating system')
end

cd(ScriptFolder)

%Compile astronomical functions that use the SOFA code
%Compile changeEpoch
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Astronomical Code/changeEpoch.c',linkCommands{:})
%Compile starCat2Obs
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./Coordinate Systems/Shared C++ Code/','-I./','./Astronomical Code/starCat2Obs.cpp','./Coordinate Systems/Shared C++ Code/getENUAxesCPP.cpp',linkCommands{:})
%Compile aberCorr
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./Coordinate Systems/Relativity/Shared C Code/','-I./','./Astronomical Code/aberCorr.c','./Coordinate Systems/Relativity/Shared C Code/relVecAddC.c',linkCommands{:})
%Compile lightDeflectCorr
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Astronomical Code/lightDeflectCorr.c',linkCommands{:})

%Compile the atmospheric code
%Compile NRLMSISE00GasTemp
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/nrlmsise-00.20151122/','-I./','./Atmospheric Models/NRLMSISE00GasTemp.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00_data.c')
%Compile NRLMSISE00Alt4Pres
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/nrlmsise-00.20151122/','-I./','./Atmospheric Models/NRLMSISE00Alt4Pres.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00.c','./3rd_Party_Libraries/nrlmsise-00.20151122/nrlmsise-00_data.c')
%Compile simpAstroRefParam
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Atmospheric Models/simpAstroRefParam.c',linkCommands{:})

%%Compile the coordinate transforms that use the SOFA code.
%Compile GCRS2ITRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2ITRS.c',linkCommands{:})
%Compile ITRS2GCRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2GCRS.c',linkCommands{:})
%Compile GCRS2TIRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2TIRS.c',linkCommands{:})
%Compile TIRS2GCRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2GCRS.c',linkCommands{:})
%Compile TEME2ITRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TEME2ITRS.c',linkCommands{:})
%Compile ITRS2TEME
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2TEME.c',linkCommands{:})
%Compile TOD2GCRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TOD2GCRS.c',linkCommands{:})
%Compile GCRS2TOD
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2TOD.c',linkCommands{:})
%Compile MOD2GCRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/MOD2GCRS.c',linkCommands{:})
%Compile GCRS2MOD
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2MOD.c',linkCommands{:})
%Compile J2000F2ICRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/J2000F2ICRS.c',linkCommands{:})
%Compile ICRS2J2000F
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2J2000F.c',linkCommands{:})
%Compile TIRS2ITRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2ITRS.c',linkCommands{:})
%Compile ITRS2TIRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ITRS2TIRS.c',linkCommands{:})
%Compile GCRS2CIRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/GCRS2CIRS.c',linkCommands{:})
%Compile CIRS2GCRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/CIRS2GCRS.c',linkCommands{:})
%Compile CIRS2TIRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/CIRS2TIRS.c',linkCommands{:})
%Compile TIRS2CIRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/TIRS2CIRS.c',linkCommands{:})
%Compile G2ICRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/G2ICRS.c',linkCommands{:})
%Compile ICRS2G
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2G.c',linkCommands{:})
%Compile ICRS2Ecliptic
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ICRS2Ecliptic.c',linkCommands{:})
%Compile ecliptic2ICRS
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Celestial and Terrestrial Systems/ecliptic2ICRS.c',linkCommands{:})

%Compile ellips2Cart
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/sofa/src/','./Coordinate Systems/ellips2Cart.c',linkCommands{:})

%Compile the time conversion functions that use the SOFA code.
%Compile UTC2Cal
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/UTC2Cal.c','Coordinate Systems/Time/Shared C Code/UTC2CalC.c',linkCommands{:})
%Compile UTC2DayCount
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/UTC2DayCount.c','Coordinate Systems/Time/Shared C Code/UTC2DayCountC.c',linkCommands{:})
%Compile Cal2UTC
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/Cal2UTC.c',linkCommands{:})
%Compile UTC2TAI
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/UTC2TAI.c',linkCommands{:})
%Compile TAI2TT
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TAI2TT.c',linkCommands{:})
%Compile TCG2TT
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TCG2TT.c',linkCommands{:})
%Compile TT2TAI
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2TAI.c',linkCommands{:})
%Compile TT2TCG
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2TCG.c',linkCommands{:})
%Compile TT2GMST
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2GMST.c',linkCommands{:})
%Compile TT2GAST
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TT2GAST.c',linkCommands{:})
%Compile TAI2UTC
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TAI2UTC.c',linkCommands{:})
%Compile BesselEpoch2TDB
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/BesselEpoch2TDB.c',linkCommands{:})
%Compile TDB2BesselEpoch
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TDB2BesselEpoch.c',linkCommands{:})
%Compile TDB2TCB
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TDB2TCB.c',linkCommands{:})
%Compile TCB2TDB
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/TCB2TDB.c',linkCommands{:})
%Compile cumLeapSec
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/cumLeapSec.c',linkCommands{:})
%Compile JulDate2JulEpoch
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/JulDate2JulEpoch.c',linkCommands{:})
%Compile JulEpoch2JulDate
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','./Coordinate Systems/Time/JulEpoch2JulDate.c',linkCommands{:})
%Compile TT2UT1
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2UT1.c',linkCommands{:})
%Compile TAI2UT1
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TAI2UT1.c',linkCommands{:})
%Compile TT2TCB
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2TCB.c',linkCommands{:})
%Compile TT2TDB
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TT2TDB.c',linkCommands{:})
%Compile TDB2TT
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/sofa/src/','-I./','-I./Mathematical Functions/Shared C Code/','-I./Coordinate Systems/Time/Shared C Code/','./Coordinate Systems/Time/TDB2TT.c',linkCommands{:})

%%Compile other astronomical code
%Compile approxSolarSysVec
mex('-v','CFLAGS="$CFLAGS -std=c99"','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./','-I./3rd_Party_Libraries/sofa/src/','./Astronomical Code/approxSolarSysVec.c',linkCommands{:})
%Compile propagateOrbitSGP4
mex('-v','-largeArrayDims','-U__STDC_UTF_16__','-outdir','./0_Compiled_Code/','-I./3rd_Party_Libraries/SGP4/cpp/','-I./','./Astronomical Code/propagateOrbitSGP4.cpp','./3rd_Party_Libraries/SGP4/cpp/sgp4unit.cpp')

%Get rid of the library as nothing else is going to link to it.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('rm ./0_Compiled_Code/libsofa_c.a');
else%Windows
    if(strcmp(compName,'gcc'))
        system('del .\0_Compiled_Code\libsofa_c.a');
    else
        system('del .\0_Compiled_Code\libsofa_c.lib');
    end
end

%Compile the AIS library
cd ./3rd_Party_Libraries/libais-0048ceb/src/libais/
%Run commands on the command line to build the library and place it in a
%known location.
if(isunix()||ismac())%*NIX/ Mac OS X
    system('make -f Makefile-custom libais.a');
    system('mv ./libais.a ../../../../0_Compiled_Code/libais.a');
    system('make -f Makefile-custom clean');
    linkCommands{1}='./0_Compiled_Code/libais.a';
elseif(ispc())%Windows    
    compPath=getenv('COMPILER');
    [~,compName]=fileparts(compPath);
    cc=[compPath,' ',getenv('COMPFLAGS'),' -std=c++11',' ','-U__STDC_UTF_16__ '];
    ll=getenv('LINKER');

    %Create object files for all of the C++ files that are used.
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
    
    if(strcmp(compName,'gcc'))%If here, then minGW is probably being used.
        %Store the names of the object files
        system('dir /b *.o > temp.lst');
        
        %Link the object files into an archive.
        system('ar rv ./libais.a @temp.lst');
        system('ranlib ./libais.a');
        
        %Move the library into the folder with the compiled code.
        system('move libais.a ..\..\..\..\0_Compiled_Code\libais.a');
        
        %Get rid of the object files
        system('del *.o');
        
        %Get rid of the list of names of the object files.
        system('del temp.lst');
         linkCommands{1}='./0_Compiled_Code/libais.a'; 
    else    
        %Store the names of the object files
        system('dir /b *.obj > temp.lst');
        %Combine the object files into a library. it is assumed that the linker
        %takes options akin to the link function in Microsoft's SDK, so that we
        %can tell it to make things into a library rather than trying to
        %compile it into a program.
        system([ll,' -lib /nologo /out:libais.lib @temp.lst']);
        %Move the library into the folder with the compiled code.
        system('move libais.lib ..\..\..\..\0_Compiled_Code\libais.lib');
        %Get rid of the object files
        system('del *.obj');
        %Get rid of the list of names of the object files.
        system('del temp.lst');
        linkCommands{1}='-llibais';
        linkCommands{2}='-L./0_Compiled_Code';
    end
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
    if(strcmp(compName,'gcc'))
        system('del .\0_Compiled_Code\libais.a');
    else
        system('del .\0_Compiled_Code\libais.lib');
    end
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
