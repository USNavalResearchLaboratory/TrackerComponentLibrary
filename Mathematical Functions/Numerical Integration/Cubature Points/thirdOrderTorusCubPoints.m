function [xi,w]=thirdOrderTorusCubPoints(rho1,rho2)
%%THIRDORDERTORUSCUBPOINTS Generate third-order cubature points for
%             integration over a 3D torus with a cirular cross section.
%             The region is specified as 
%             (sqrt(x^2+y^2)-rho1)^2+z^2<=rho2^2 with 0<rho2<=rho1<Inf.
%
%INPUTS: rho1,rho2 The two defining parameters of the torus as given above
%              with 0<rho2<=rho1<Inf. rho2 is the radius of the filled-in
%              portion of the ring and rho1 is the distance of the center
%              of the ring (the hole) from the center of the filled-in
%              portion of the ring.
%
%OUTPUTS: xi A 3X8 matrix containing the cubature
%            points (Each "point" is a vector).
%          w A 8X1 vector of the weights associated with
%            the cubature points.
%
%This function implements formula Tor_{3}:S_2: 3-1 on page 340 of [1], 8
%points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

r=sqrt(rho1^2+(3/4)*rho2^2);
s=rho2/2;

xi=[PMCombos([r;0;s]),PMCombos([0;r;s])];

V=2*pi^2*rho1*rho2^2;
B=(1/8)*V;

w=B*ones(8,1);

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
