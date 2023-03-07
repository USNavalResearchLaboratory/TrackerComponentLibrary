function Area=unitVectors2SpherTriangArea(u,r)
%%UNITVECTORS2SPHERTRIANGAREA Given unit vectors pointing from the center
%     of a 3D sphere to three points on the sphere, as well as the radius
%     of the sphere, determine the surface area of a spherical triangle.
%
%INPUTS: u A 3X3 matrix such that u(:,i) is the ith unit vector pointing
%          from the center of the sphere to the surface of the sphere.
%        r The scalar radius of the sphere. If this is omitted or an empty
%          matrix is passed, the default of 1 is used.
%
%OUTPUTS: Area The area of the specified triangle on the surface of a
%              sphere of radius r.
%
%Equation 1 in [1] is implemented. The absolute value of the atan function
%is taken so that the result is always positive.
%
%EXAMPLE:
%We find the surface area of a quadrant of a unit sphere. Since the surface
%area of an entire sphere is 4*pi*r^2, the area of a quadrant will be pi/2.
%The relative error is zero.
% u=[[1;0;0],[0;1;0],[0;0;1]];
% r=1;
% Area=unitVectors2SpherTriangArea(u,r);
% trueArea=pi/2;
% RelErr=(Area-trueArea)./trueArea
%
%REFERENCES:
% [1] F. Eriksson, "On the measure of solid angles," Mathematics Magazine,
%     vol. 63, no. 3, pp. 184-187, 1990.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(r))
    r=1;
end

a=u(:,1);
b=u(:,2);
c=u(:,3);

Area=2*r^2*abs(atan(dot(a,cross(b,c))/(1+dot(b,c)+dot(c,a)+dot(a,b))));

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
