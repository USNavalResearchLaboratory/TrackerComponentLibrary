function [xi,w]=firstOrderNDimCubPoints(numDim,algorithm)
%%FIRSTORDERNDIMCUBPOINTS Generate first-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>0.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are:
%               0 Formula Cn 1-1 in [1], pg. 229, 1 point.
%               1 (the default if omtted or an empty matrix is passed)
%                 Formula Cn 1-2 in [1], pg. 229. 2^numDim points.
%               2 Formula C2 1-1 in [1], pg. 242, 4 points, numDim=2.
%               3 Formula C2 1-2 in [1], pg. 243, 9 points, numDim=2.
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=1; 
end

switch(algorithm)
    case 0%Cn 1-1, 1 point in [1].
        V=2^numDim;
        xi=zeros(numDim,1);
        w=V;
    case 1%Cn 1-2 2^n points in [1].
        V=2^numDim;
        xi=PMCombos(ones(numDim,1));
        w=1/2^numDim*ones(numDim,1)*V;
    case 2%C2 1-1 in [1], pg. 242, 4 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        V=2^numDim;
        w=(1/4)*V*ones(4,1);
        xi=PMCombos([1;1]);
    case 3%C2 1-2 in [1], pg. 243, 9 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        V=2^numDim;
        w=[250/225*V;-8/225*V*ones(4,1);7/900*V*ones(4,1)];
        xi=[[0;0],fullSymPerms([1;0]),PMCombos([1;1])];
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
