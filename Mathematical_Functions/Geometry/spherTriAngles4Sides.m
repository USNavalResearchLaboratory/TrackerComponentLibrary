function [A,B,C]=spherTriAngles4Sides(a,b,c,r)
%%SPHERTRIANGLES4SIDES Given the lengths of the sides of a spherical
%       triangle, find the angles of the triangle. Note that the angles of
%       a spherical triangle do not sum to pi. This is just an application
%       of the spherical law of cosines. Note that unlike with a EUclidean
%       triangle, one cannot just enter arbitrary distances for a, b, and c
%       (even if much smaller than the circumference of the sphere),
%       because certain distance relations are not possible and will result
%       in non-real solutions.
%
%INPUTS: a, b, c The lengths of the sides of the spherical triangle
%                in any order. If one wishes to convert multiple triangles
%                at once, these can be matrices.
%              r The scalar radius of the sphere on which the distances are
%                measured. If omitted or an empty matrix is passed, the
%                default of 1 is used.
%
%OUTPUTS: A, B, C The angles of the spherical triangle in radians. A is the
%                 angle opposite side a, B the angle opposite side b and C
%                 the angle opposite side C. The dimensions correspond to
%                 those of the inputs.
%
%The spherical law of cosines for the sides of a spherical triangle is
%stated in [1].
%
%EXAMPLE:
%First, we choose 3 points on the sphere, find the distances between them
%and then show that we can get the angles. Then, we just set some arbitrary
%distances that do not form a possible triangle on the spehre and we see
%that the angular results are complex.
%Choose 3 points in term of latitude-longitude on a sphere.
% llPoints=deg2rad([19,  21, 27;
%                 -158,-147,-152]);
% r=1;
% a=greatCircleDistance(llPoints(:,1),llPoints(:,2),r);
% b=greatCircleDistance(llPoints(:,2),llPoints(:,3),r);
% c=greatCircleDistance(llPoints(:,3),llPoints(:,1),r);
% [A,B,C]=spherTriAngles4Sides(a,b,c,r);
% ADeg=rad2deg(A)
% BDeg=rad2deg(B)
% CDeg=rad2deg(C)
% a=1;
% b=2;
% c=3;
% [A,B,C]=spherTriAngles4Sides(a,b,c,r)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Law of Cosines." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/LawofCosines.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(r))
    r=1;
end

a=a./r;
b=b./r;
c=c./r;

cosa=cos(a);
sina=sin(a);
cosb=cos(b);
sinb=sin(b);
cosc=cos(c);
sinc=sin(c);

A=acos((cosa-cosb.*cosc)./(sinb.*sinc));
B=acos((cosb-cosa.*cosc)./(sina.*sinc));
C=acos((cosc-cosa.*cosb)./(sina.*sinb));

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
