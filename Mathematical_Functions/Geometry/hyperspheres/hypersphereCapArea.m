function area=hypersphereCapArea(r,numDim,lat)
%%HYPERSPHERECAPAREA Determine the surface area of a sphere in a given
%               number of dimensions that has been cut by a plane leaving
%               only a cap. This function returns the surface area of the
%               cap (the smaller half of the cut sphere), only considering
%               the surface on the outside of the hypersphere.
%
%INPUTS:     r The radius of the hypersphere. r>0.
%       numDim The number of dimensions of the hypersphere. numDim>=1. One
%              dimensions is a line, two a circle, three a sphere, etc.
%          lat Defining the cap as beginning some angle above the equator
%              of the sphere, this is the angle above the equator in
%              radians 0<=lat<=pi/2.
%
%OUTPUTS: area The surface area of the hyperspherical cap.
%
%The formula is given in [1]. For an illustration of a hyperspherical cap,
%see the image in [2].
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

n=numDim;
phi=pi/2-lat;
area=(1/2)*hypersphereArea(r,n)*betainc(sin(phi)^2,(n-1)/2,1/2);

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
