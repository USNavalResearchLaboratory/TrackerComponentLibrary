function Area=sideAngles2SpherTriangArea(sideAng,r)
%%SIDEANGLES2SPHERTRIANGAREA Given the length of sides of a triangle on a
%          sphere (scaled to a unit sphere) -these are angles between unit
%          vectors from the center to the points on any sphere, determine
%          the surface area of the triangle on a sphere of radius r. This
%          is similar to the spherTriangArea4Sides function.
%
%INPUTS: sideAngs The distances across the surface of a unit sphere between
%          each of the three points in the triangle. Such distances can be
%          computed using greatCircleDistance. These are angles between
%          unit vectors from the center of the sphere to each of the points
%          om the triangle.
%        r The scalar radius of the sphere. If this is omitted or an empty
%          matrix is passed, the default of 1 is used.
%
%OUTPUTS: Area The area of the specified triangle on the surface of a
%              sphere of radius r.
%
%Equations 2 and 3 in [1] are used. 
%
%EXAMPLE:
%Here, we find the great circle distance between 3 points on the sphere and
%compute the area of the spherical triangle. The result is compared to the
%output of unitVectors2SpherTriangArea and the output of this function and
%that function are shown to be the same, within finite precision limits.
% r=20;%Sphere radius
% llVals=deg2rad([19.129599,  20.452749,  39.491324;
%                -155.251734,-139.069286,-146.808718]);
% %a=1,f=0 means unit vectors on a sphere.
% unitVectors=ellips2Cart(llVals,1,0);
% sideAngs=zeros(3,1);
% sideAngs(1)=greatCircleDistance(llVals(:,1),llVals(:,2),1);
% sideAngs(2)=greatCircleDistance(llVals(:,2),llVals(:,3),1);
% sideAngs(3)=greatCircleDistance(llVals(:,3),llVals(:,1),1);
% Area1=sideAngles2SpherTriangArea(sideAngs,r)
% Area2=unitVectors2SpherTriangArea(unitVectors,r)
% RelDiff=(Area1-Area2)/Area2
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

cosA=cos(sideAng(1));
cosB=cos(sideAng(2));
cosC=cos(sideAng(3));

%Equation 3.
P=sqrt(1-cosA^2-cosB^2-cosC^2+2*cosA*cosB*cosC);

%Equation 2.
Area=2*r^2*atan(P/(1+cosA+cosB+cosC));

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
