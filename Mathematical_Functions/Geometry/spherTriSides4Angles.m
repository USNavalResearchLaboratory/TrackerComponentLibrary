function [a,b,c]=spherTriSides4Angles(A,B,C,r)
%%SPHERTRISIDES4ANGLES Given the angles of a spherical triangle and the
%               radius of the sphere on which the triangle is inscribed,
%               determine the lengths of the sides of the triangle. 
%
%INPUTS: A, B, C The angles of the spherical triangle in radians in any
%                order. If one wishes to convert multiple triangles at
%                once, these can be matrices.
%
%INPUTS: a, b, c The scalar lengths of the sides of the spherical triangle.
%                A is the angle opposite side a, B the angle opposite side
%                b and C the angle opposite side C. The dimensions are
%                consistent with the dimensions of the inputs.
%
%The spherical law of cosines for the angles of a spherical triangle is
%stated in [1].
%
%EXAMPLE:
%Given three angles, we find the sides of a triangle using this function.
%Then, we show that the result is correct, because we can invert the
%solution using the spherTriAngles4Sides function. The relative error is
%shown and is zero within finite precision limits.
% %The angular sum must be >180 degrees (pi radians). The spherical
% %triangle area is proportional to the sum of the angles.
% A=deg2rad(35);
% B=deg2rad(65);
% C=deg2rad(110);
% r=2;
% [a,b,c]=spherTriSides4Angles(A,B,C,r);
% [ABack,BBack,CBack]=spherTriAngles4Sides(a,b,c,r);
% RelErr=max(abs(([ABack,BBack,CBack]-[A,B,C])./[A,B,C]))
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

cosA=cos(A);
sinA=sin(A);
cosB=cos(B);
sinB=sin(B);
cosC=cos(C);
sinC=sin(C);

a=r.*acos((cosA+cosB.*cosC)./(sinB.*sinC));
b=r.*acos((cosB+cosA.*cosC)./(sinA.*sinC));
c=r.*acos((cosA.*cosB+cosC)./(sinA.*sinB));

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
