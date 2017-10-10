function zC=bipolar2Cart(zBi,a)
%%BIPOLAR2CART Convert from the 2D bipolar coordinate system that is based
%              on Apollonian circles to 2D Cartesian coordinates. Bipolar
%              coordinates are often used when solving certain artial
%              differential equations (such as Laplace's equation in 2D)
%              where the coordinate system allows a separation of
%              variables.
%
%INPUTS: zBi A 2XN set of bipolar [nu;u] points to convert into Cartesian
%            coordinates. The nu values map to unique Cartesian locations
%            for all real numbers. The u values map periodically to
%            Cartesian locations, periodic with period 2*pi.
%          a An arbitrary positive constant that defines the coordinate
%            system. If this parameter is omitted or an empty matrix is
%            passed, then the default of 1 is used.
%
%OUTPUTS: zC The 2XN set of values in zBi converted intp [x;y] Cartesian
%            points.
%
%The coordinate system is described in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Bipolar Coordinates." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/BipolarCoordinates.html
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(a))
    a=1; 
end

N=size(zBi,2);

zC=zeros(2,N);

nu=zBi(1,:);
u=zBi(2,:);

denom=cosh(nu)-cos(u);

zC(1,:)=a*sinh(nu)./denom;%x
zC(2,:)=a*sin(u)./denom;%y

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
