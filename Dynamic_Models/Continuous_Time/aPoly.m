function [aVal,aJacob,aHess,papt]=aPoly(x,numDim)
%%APOLY The drift function for a linear continuous-time motion model of a
%       given order in a specified number of Cartesian dimensions. The
%       order of the linear filter, that is the number of moments of
%       position, does not need to be explicitly specified.
%
%INPUTS: x The xDimXN state vector of N targets in the order of
%          [position;velocity;acceleration;etc] for however many
%          derivatives of position there are.
%   numDim The number of dimensions of the simulation problem. If the
%          numDim parameter is omitted, then numDim=3 (3D motion) is
%          assumed. The dimensionality of the state must be an integer
%          multiple of numDim.
%
%OUTPUTS: aVal The xDimXN time-derivative of the N state vectors under the
%              linear motion model.
%       aJacob The xDimXxDim Jacobian of aVal (it is the same for all x and
%              is not repeated N times). This is such that aJacob(:,k) is
%              the partial derivative of aVal with respect to the kth
%              element of x.
%        aHess The xDimXxDimXxDim hypermatrix such that aHess(:,k1,k2) is
%              the second partial derivative of aVal with respect to
%              elements k1 and k2 of x (all zero in this instance). It is
%              the same for all x.
%         papt The xDimX1 partial derivative of aVal with respect to time
%              (all zero in this instance). It is the same for all x.
%
%The drift function corresponds to the state transition given in
%discrete-time by the function FPolyKal.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    numDim=3; 
end

numTar=size(x,2);
xDim=size(x,1);

aVal=[x((numDim+1):(xDim),:);zeros(numDim,numTar)];

if(nargout>1)
    aJacob=[zeros(xDim,numDim),[eye(xDim-numDim,xDim-numDim);zeros(numDim,xDim-numDim)],zeros(xDim,numDim)];

    if(nargout>2)
        aHess=zeros(xDim,xDim,xDim);

        if(nargout>3)
            papt=zeros(xDim,1);
        end
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
