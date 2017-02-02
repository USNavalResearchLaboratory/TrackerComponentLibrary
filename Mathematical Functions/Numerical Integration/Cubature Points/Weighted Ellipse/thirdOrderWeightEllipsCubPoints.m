function [xi,w]=thirdOrderWeightEllipsCubPoints(a,b,algorithm)
%%THIRDORDERWEIGHTELLIPSCUBPOINTS Generating third-order cubature points
%               for integrating over a 2D elliptical region (aligned with
%               the coordinate axes, x and y) having weighting function
%               w(x)=1/(sqrt((x-c)^2+y^2)*sqrt((x+c)^2+y^2) where
%               c=sqrt(a^2-b^2) and a and b are the parameters of the
%               ellipse such that x^2/a^2+y^2/b^2=1 with b<=a.
%
%INPUTS: a,b The positive, real, scalar parameters of the ellipse such that
%            x^2/a^2+y^2/b^2=1 defines the ellipse and b<=a.
%  algorithm A parameter selecting the algorithm to use for the points.
%            Possible values are:
%            0 (the default if omitted or an empty matrix is passed)
%              Formula ELP 3-1 in [1], pg. 336, 4 points.
%            1 Formula ELP 3-2 in [1], pg. 336, 4 points.
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
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(b>a)
    error('b must be <=a')
end

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

c=sqrt(a^2-b^2);
V=2*pi*log((a+b)/c);

switch(algorithm)
    case 0%ELP 3-1 in [1], pg. 336, 4 points.
        r=sqrt((pi/V)*(a*b+(1/(2*pi))*c^2*V));
        s=sqrt((pi/V)*(a*b-(1/(2*pi))*c^2*V));
        
        B=V/4;
        
        xi=[PMCombos([r;0]),PMCombos([0;s])];
        w=B*ones(4,1);
    case 1%ELP 3-2 in [1], pg. 336, 4 points.
        r=sqrt((pi/(2*V))*(a*b+(1/(2*pi))*c^2*V));
        s=sqrt((pi/(2*V))*(a*b-(1/(2*pi))*c^2*V));
        
        xi=PMCombos([r;s]);
        w=(V/4)*ones(4,1);
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
