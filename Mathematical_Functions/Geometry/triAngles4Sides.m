function [A,B,C]=triAngles4Sides(a,b,c)
%%TRIANGLES4SIDES Given the sides of a triangle, find the angles of the
%         triangle. This just applies the law of cosines.
%
%INPUTS: a, b, c The scalar lengths of the sides of the triangle in any
%                order.
%
%OUTPUTS: A, B, C The angles of the triangle in radians. A is the angle
%                 opposite side a, B the angle opposite side b and C the
%                 angle opposite side C.
%
%The law of cosines is given in [1], among many other places.
%
%EXAMPLE:
%This recovers the angles of a 30-60-90 traingle.
% a=1;
% b=sqrt(3);
% c=2;
% [A,B,C]=triAngles4Sides(a,b,c);
% ADeg=rad2deg(A)
% BDeg=rad2deg(B)
% CDeg=rad2deg(C)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Law of Cosines." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/LawofCosines.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

a2=a.*a;
b2=b.*b;
c2=c.*c;

A=acos((-a2+b2+c2)./(2*b.*c));
B=acos((a2-b2+c2)./(2*a.*c));
C=pi-A-B;
%We could have used:
%C=acos((a2+b2-c2)./(2*a.*b));
%But we know that the angles in a triangle in Euclidean geometry sum to 180
%degrees.

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
