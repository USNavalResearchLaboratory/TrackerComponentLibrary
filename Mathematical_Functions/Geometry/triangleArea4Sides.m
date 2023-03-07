function A=triangleArea4Sides(a,b,c)
%%TRIANGLEAREA4SIDES Find the area of a triangle given the lengths of the 3
%                    sides. Heron's formula as in [1] is used.
%
%INPUTS: a,b,c The positive lengths of the 3 sides of the triangle. These
%              can be matrices if one wishes to evaluate multiple triangles
%              at once.
%
%OUTPUTS: A The area of the triangle(s). The dimensions of A match those of
%           the inputs.
%
%EXAMPLE:
%Here, we have a 24-30-18 triangle. Its area is 216.
% A=triangleArea4Sides(24,30,18)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Heron's Formula." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/HeronsFormula.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    s=(a+b+c)./2;
    A=sqrt(s.*(s-a).*(s-b).*(s-c));

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
