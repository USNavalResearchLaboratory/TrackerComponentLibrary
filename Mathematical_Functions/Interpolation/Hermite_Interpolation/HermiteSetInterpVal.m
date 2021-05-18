function yInterp=HermiteSetInterpVal(pInterp,a,c,xRange,numDerivs)
%%HERMITESETINTERPVAL Given a set of interpolating polynomial coefficients
%           a and control points c for Hermite interpolation from the
%           function findHermiteInterpPolySet, this function evaluates the
%           interpolated values at the points in pInterp.
%
%INPUTS: pInterp A numPointsX1 or 1XnumPoints set of points at which the
%           interpolated values are desired.
%         a The numDimsXnumCoeffXnumInterpRegions
%           hypermatrix of the interpolating polynomial coefficients for
%           each dimension, moment, and region. This and the following two
%           other inputs can be obtained using the findHermiteInterpPolySet
%           function.
%         c The (numDims*numMoments2Interp)X(numCoeff-1)XnumInterpRegions
%           hypermatrix of control points for the interpolating
%           polynomials.
%    xRange The length (numInterpRegions+1) vector of boundary points of
%           the interpolating regions.
% numDerivs The number of derivatives to interpolate. A value of zero means
%           only interpolate position. 1 means velocity, 2 acceleration,
%           etc.
%
%OUTPUTS: yInterp The (numDims*(numDerivs+1))XnumPoints set of interpolated
%                 state values at the given points. All dimensions for a
%                 particular deirvative value are given before later ones.
%                 for example, for position and velocity in 3D, the
%                 ordering is [x;y;z;xDot;yDot;zDot].
%
%This function uses binSearch on each point to find the region in xRange
%indicating which values in a and c to use for the interpolation. The
%function then performs interpolation using multiple calls to
%polyValNewton. See the sample code in the comments to
%findHermiteInterpPolySet for an example of how this function is used.
%
%May 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numDerivs))
    numDerivs=0; 
end

numDims=size(a,1);
numPoints=length(pInterp);
numRanges=length(xRange)-1;

yInterp=zeros(numDims*(numDerivs+1),numPoints);

if(numDerivs==0)
    for curPoint=1:numPoints
        %Determine the start of the relevant interpolating range.
        [~,idx]=binSearch(xRange,pInterp(curPoint),1);
        idx=min(idx,numRanges);

        yInterp(:,curPoint)=polyValNewton(pInterp(curPoint),a(:,:,idx),c(:,:,idx));
    end
else
    for curPoint=1:numPoints
        %Determine the start of the relevant interpolating range.
        [~,idx]=binSearch(xRange,pInterp(curPoint),1);
        idx=min(idx,numRanges);
        
        [pDer,p]=polyDerValNewton(pInterp(curPoint),a(:,:,idx),c(:,:,idx),numDerivs);
        yInterp(:,curPoint)=[p;reshape(pDer,[numDims*numDerivs,1])];
    end
end
    
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
