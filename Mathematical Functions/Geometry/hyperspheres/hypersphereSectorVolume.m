function vol=hypersphereSectorVolume(r,numDim,lat)
%%HYPERSPHERECAPVOLUME Determine the volume of a sphere in a given number
%               of dimensions that has been cut by a particular solid angle
%               leaving a cone and a cap (a spherical sector).
%
%INPUTS:     r The radius of the hypersphere. r>0.
%       numDim The number of dimensions of the hypersphere. numDim>=1. One
%              dimensions is a line, two a circle, three a sphere, etc.
%          lat The angle with which the cone starts above the equator in
%              radians 0<=lat<=pi/2.
%
%OUTPUTS: vol The volume of the specified hyperspherical sector.
%
%The formula is given in [1]. For an illustration of a hyperspherical
%sector, see the image in [2].
%
%EXAMPLE:
%An expression for the volume of a spherical sector in 3D is given as
%vol=2/3*pi*r^2*h;  where h=r-r*sin(lat);
%Thus, we can try it out:
% r=5;
% numDim=3;
% lat=0.4;
% h=r-r*sin(lat);
% vol1=2/3*pi*r^2*h
% vol2=hypersphereSectorVolume(r,numDim,lat)
%One finds that vol1 and vol2 are the same.
%
%REFERENCES:
%[1] S. Li, "Concise formulas for the area and volume of a hyperspherical
%    cap," Asian Journal of Mathematics and Statistics, vol. 4, no. 1, pp.
%    66-70, 2011.
%[2] Weisstein, Eric W. "Spherical Sector." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/SphericalSector.html
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=numDim;
phi=pi/2-lat;

vol=(1/2)*hypersphereVolume(r,n)*betainc(sin(phi)^2,(n-1)/2,1/2);

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