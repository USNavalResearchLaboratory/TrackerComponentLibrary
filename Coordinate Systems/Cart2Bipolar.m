function zBi=Cart2Bipolar(zC,a)
%%CART2BIPOLAR Convert from 2D Cartesian coordinates into the bipolar
%              coordinate system that is based on Apollonian circles.
%              Bipolar coordinates are often used when solving certain
%              partial differential equations (such as Laplace's equation
%              in 2D) where the coordinate system allows a separation of
%              variables.
%
%INPUTS: zC A 2XN set of Cartesian [x;y] points to convert into bipolar
%           coordinates.
%         a An arbitrary positive constant that defines the coordinate
%           system. If this parameter is omitted or an empty matrix is
%           passed, then the default of 1 is used.
%
%OUTPUTS: zBi The 2XN set of values in aZ converted into bipolar
%             coordinates [nu;u]. The nu values range from -Inf to Inf and
%             the u values range from -pi to pi.
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

N=size(zC,2);

zBi=zeros(2,N);

x=zC(1,:);
y=zC(2,:);

zBi(1,:)=atanh(2*a*x./(a^2+x.^2+y.^2));%nu
zBi(2,:)=atan2(2*a*y,x.^2+y.^2-a^2);%u

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
