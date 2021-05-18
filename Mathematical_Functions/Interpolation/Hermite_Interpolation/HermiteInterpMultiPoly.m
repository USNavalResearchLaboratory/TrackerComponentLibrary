function [a,c]=HermiteInterpMultiPoly(x,y,numDims)
%%HERMITEINTERPMULTIPOLYS Given a multidimensional function y(x) evaluated
%             a certain values of x (control points) along with an
%             arbitrary number of derivatives of y(x) at those points, this
%             function returns the coefficients of the Hermite
%             interpolating polynomials fitting those points across all
%             dimensions. This function differs from
%             findHermiteInterpPolySet in that this function tries to find
%             one high-order Hermite interpolation polynomial for each
%             dimension that matches all of the points given, whereas
%             findHermiteInterpPolySet finds a set of low-order
%             polynomials, which is more practical for large regions as
%             finite precision errors can dominate this function. Compare
%             this to the function HermiteInterpPoly, which can only handle
%             scalar values.
%
%INPUTS: x An NpX1 or 1XNp vector of real, scalar values at which the
%          function y and its derivative are given. The values are assumed
%          to be given in increasing order. Np>=2.
%        y A (numMoments*numDims)XNp matrix of values of the
%          multidimensional function y(x) and its derivatives evaluated at
%          the points in x. All values for a particular derivative are
%          given before any of the next derivative. For example, to
%          interpolate in 3D position and velocity, the ordering would be
%          [x;y;z;xDot;yDot;zDot].
%  numDims The scalar number of dimensions present. This will typically be
%          2 or 3 when dealing with target states. If this parameter is
%          omitted or an empty matrix is passed, then numDims=1 is used.
%
%OUTPUTS: a The numDimsXnumCoeff matrix of the interpolating polynomial
%           coefficients for each dimension and moment. This and the
%           next output can be used in the function polyValNewton to
%           interpolate values.
%         c The numDimsX(numCoeff-1) matrix of control points for the
%           interpolating polynomials.
%
%This function calls HermiteInterpPoly for each of the dimensions present
%and stacks the results by dimension.
%
%May 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(numDims))
    numDims=1;
end

if(numDims==1)
    [a,c]=HermiteInterpPoly(x,y);
    return
end

Np=length(x);
stateSize=size(y,1);

numMatch=stateSize/numDims;
numCoeff=numMatch*Np;

a=zeros(numDims,numCoeff);
c=zeros(numDims,numCoeff-1);

for curDim=1:numDims
    yCur=y(curDim:numDims:stateSize,:);
    [a(curDim,:),c(curDim,:)]=HermiteInterpPoly(x,yCur);
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
