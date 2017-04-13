function vol=hypersphereCapVolume(r,numDim,lat)
%%HYPERSPHERECAPVOLUME Determine the volume of a sphere in a given number
%               of dimensions that has been cut by a plane leaving only a
%               cap. This function returns the volume of the smaller half.
%
%INPUTS:     r The radius of the hypersphere. r>0.
%       numDim The number of dimensions of the hypersphere. numDim>=1. One
%              dimensions is a line, two a circle, three a sphere, etc.
%          lat Defining the cap as beginning some angle above the equator
%              of the sphere, this is the angle above the equator in
%              radians 0<=lat<=pi/2.
%
%OUTPUTS: volume The volume of the hyperspherical cap.
%
%The formula is given in [1]. For an illustration of a hyperspherical cap,
%see the image in [2].
%
%EXAMPLE:
%In three dimensions, the volume of a hyperspherical cap is
%vol=1/3*pi*r^3*(2-3*sin(lat)+sin(lat)^3) as in [2].
%Comparing it to this function
% lat=0.4;
% r=3;
% vol1=1/3*pi*r^3*(2-3*sin(lat)+sin(lat)^3)
% vol2=hypersphereCapVolume(r,3,lat)
%One obtains equivalent solutions.
%
%REFERENCES:
%[1] S. Li, "Concise formulas for the area and volume of a hyperspherical
%    cap," Asian Journal of Mathematics and Statistics, vol. 4, no. 1, pp.
%    66-70, 2011.
%[2] Weisstein, Eric W. "Spherical Cap." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/SphericalCap.html
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

phi=pi/2-lat;

n=numDim;
vol=(1/2)*hypersphereVolume(r,n)*betainc(sin(phi)^2,(n+1)/2,1/2);
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
