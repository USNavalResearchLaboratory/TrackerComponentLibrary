function [xi,w]=thirdOrderHexagonCubPoints(algorithm)
%%THIRDORDERHEXAGONCUBPOINTS Generate third-order cubature points for
%               integration over a two-dimensional hexagon with vertices
%               (+/-1,0) and (+/-(1/2),+/-(1/2)*sqrt(3)).
%
%INPUTS: algorithm An optional parameter that selects the algorithm to use.
%                  Possible values are:
%                  0 (The default if omited or an empty matrix is passed)
%                    Formula H_2 3-1 in [1], pg. 335, 4 points.
%                  1 H_2 3-2 in [1], pg. 336, 7 points.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%H_2 3-1 in [1], pg. 335, 4 points.
        r=sqrt((5/12));
        xi=fullSymPerms([r;0]);
        
        V=3*sqrt(3)/2;
        w=(1/4)*V*ones(4,1);
    case 1%H_2 3-2 in [1], pg. 336, 7 points.
        r=1;
        s=(1/2);
        t=sqrt(3/4);
        
        xi=[[0;0],PMCombos([r;0]),PMCombos([s;t])];
        
        V=3*sqrt(3)/2;
        A=(42/72)*V;
        B=(5/72)*V;
        
        w=[A;B*ones(6,1)];
    otherwise
        error('Unknown algorithm specified')
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
