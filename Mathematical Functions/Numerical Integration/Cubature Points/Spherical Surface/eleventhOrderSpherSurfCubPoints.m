function [xi,w]=eleventhOrderSpherSurfCubPoints(numDim,algorithm)
%%ELEVENTHORDERSPHERSURFCUBPOINTS Generate eleventh-order cubature
%               points for integration over the surface of a unit (hyper-)
%               sphere (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points
%               to be generated. Currently, only numDim=3 is supported.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula U3 11-1 in [1], pg. 301, 50 points, numDim=3.
%               1 Formula U3 11-2 in [1], pg. 301, 56 points, numDim=3.
%               2 Formula U3 11-3 in [1], pg. 301, 62 points, numDim=3.
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

if(numDim~=3)
   error('Only 3D points are supported'); 
end

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%U3 11-1 in [1], pg. 301, 50 points, numDim=3.
        V=2*pi^(numDim/2)/gamma(numDim/2);
        r=1;
        s=1/sqrt(2);
        t=1/sqrt(3);
        u=1/sqrt(11);
        v=sqrt(9/11);
        B1=9216*V/725760;
        B2=16384*V/725760;
        B3=15309*V/725760;
        B4=14641*V/725760;
        
        xi=[fullSymPerms([r;0;0]),fullSymPerms([s;s;0]),PMCombos([t;t;t]),fullSymPerms([u;u;v])];
        w=[B1*ones(6,1);B2*ones(12,1);B3*ones(8,1);B4*ones(24,1)];
    case 1%U3 11-2 in [1], pg. 301, 56 points, numDim=3.
        V=2*pi^(numDim/2)/gamma(numDim/2);
        r=sqrt((9-4*sqrt(3))/33);
        s=sqrt((15+8*sqrt(3))/33);
        v=sqrt((9+4*sqrt(3))/33);
        u=sqrt((15-8*sqrt(3))/33);
        t=1/sqrt(3);
        A=V*9/560;
        B=V*(122+9*sqrt(3))/6720;
        C=V*(122-9*sqrt(3))/6720;
        
        xi=[PMCombos([t;t;t]),fullSymPerms([s;r;r]),fullSymPerms([u;v;v])];
        w=[A*ones(8,1);B*ones(24,1);C*ones(24,1)];
    case 2%U3 11-3 in [1], pg. 301, 62 points, numDim=3.
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        p=sqrt((5+sqrt(5))/10);
        q=sqrt((5-sqrt(5))/10);
        r=sqrt((3-sqrt(5))/6);
        s=sqrt((3+sqrt(5))/6);
        t=1/sqrt(3);
        u=1/2;
        v=(sqrt(5)+1)/4;
        w=(sqrt(5)-1)/4;
        B=V*625/27720;
        C=V*243/27720;
        D=V*512/27720;
        
        xi=[PMCombos([p;q;0]),PMCombos([0;p;q]),PMCombos([q;0;p]),PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([t;t;t]),fullSymPerms([1;0;0]),PMCombos([u;v;w]),PMCombos([w;u;v]),PMCombos([v;w;u])];
        w=[B*ones(12,1);C*ones(20,1);D*ones(30,1)];
    otherwise
        error('Invalid algorithm specified')
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
