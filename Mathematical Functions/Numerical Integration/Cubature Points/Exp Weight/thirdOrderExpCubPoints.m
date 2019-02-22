function [xi,w]=thirdOrderExpCubPoints(numDim,algorithm)
%%THIRDORDEREXPCUBPOINTS Generate third-order cubature points for
%               integration over real space involving the weighting
%               function w(x)=exp(-sqrt(sum(x.*x))).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%     algorithm A value indicating which algorithm should be used. Possible
%               values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula E_n^r 3-1 in [1], pg. 329, 2*numDim points.
%               1 Formula E_n^r 3-2 in [1], pg. 329, 2^numDim points.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%E_n^r 3-1 in [1], pg. 329, 2*numDim points.
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        r=sqrt(numDim*(numDim+1));
        
        w=V/(2*numDim)*ones(2*numDim,1);
        xi=fullSymPerms([r;zeros(numDim-1,1)]);
    case 1%E_n^r 3-2 in [1], pg. 329, 2^numDim points.
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        r=sqrt(numDim+1);
        
        xi=PMCombos(r*ones(numDim,1));
        w=2^(-numDim)*V*ones(2^numDim,1);
    otherwise
        error('Unknown algorithm specified');    
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
