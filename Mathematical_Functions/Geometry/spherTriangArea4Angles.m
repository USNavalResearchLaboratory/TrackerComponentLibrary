function surfArea=spherTriangArea4Angles(A,B,C,r)
%%SPHERTRIANGAREA4ANGLES Given the 3 angles of a spherical triangle as well
%                        as the radius of the sphere, find the surface area
%                        of the triangle.
%
%INPUTS: A, B, C The angles of the spherical triangle in radians in any
%                order. These can be matrices if one wishes to evaluate
%                the areas of multiple triangles at once. Spherical
%                triangles with nonzero area have a sum of angles >pi
%                radians.
%              r The radius of the sphere on which the distances are
%                measured. If omitted or an empty matrix is passed, the
%                default of 1 is used.
%
%OUTPUTS: surfArea The surface area of the spherical triangle. The
%                  dimensions of this correspond to those of the inputs.
%
%REFERENCES:
%Using equilateral triangles, we show that the area of the triangle is
%nonzero when the angles sum to opver pi and is 0 when the angles sum to pi.
%Note that a sum less than pi cannot form a valid spherical triangle.
% E=1.1*pi/3;
% r=1;
% surfArea=spherTriangArea4Angles(E,E,E,r)
% E=pi/3;
% surfArea0=spherTriangArea4Angles(E,E,E,r)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Spherical Triangle." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/SphericalTriangle.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(r))
    r=1;
end

surfArea=r^2*(A+B+C-pi);

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
