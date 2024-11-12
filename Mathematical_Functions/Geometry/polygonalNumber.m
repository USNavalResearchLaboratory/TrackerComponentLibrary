function val=polygonalNumber(n,numSides)
%%POLYGONALNUMBER Generate the nth polygonal number for a shape with
%       numSides sides. These are a type of figurate number. Polygonal
%       numbers represent the number of dots that would be needed to form a
%       shape of a filled-in regular polygon with a side length of n. THis
%       also means that one can connect the dots to form n-1 nested shapes
%       of numSides sides.
%
%INPUTS: n The length of the polygon side (=the number of nested
%          polygons+1).
% numSides The number of sides of the polygons. numSides>=3. If omitted or
%          an empty matrix is used, the default if 3 (traingular numbers).
%
%OUTPUTS: val The specified polygonal number.
%
%As an example, if n=4 and numSides=3,, then one can form a triangle on a
%grid as
% ...x...
% ..x.x..
% .x.x.x.
% x.x.x.x
%Where the xs are the points. So, as expected, polygonalNumber(4) is 10.
%For a rectangular number with a side length of 4, however, one ends up
%with a grid of
% xxxx
% xxxx
% xxxx
% xxxx
%which also corresponds to polygonalNumber(4,4) being 16.
%
%Formulae for polygonal numbers are given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Polygonal Number." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/PolygonalNumber.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==1||isempty(numSides))
    %Triangular numbers.
    numSides=3;
end

val=n*(n*(numSides-2)-numSides+4)/2;
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
