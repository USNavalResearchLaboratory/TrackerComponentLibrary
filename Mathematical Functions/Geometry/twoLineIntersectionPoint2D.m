function point=twoLineIntersectionPoint2D(line1,line2)
%%TWOLINEINTERSECTIONPOINT2D Given two (infinite) 2D lines that are not
%                     parallel, find the point of intersection. A line is
%                     specified by providing two points on the line.
%
%INPUTS:   line1 A 2X2 matrix of two points in the first line where
%                line1(1,:) are the x-coordinates and line1(2,:) are the
%                y-coordinates.
%          line2 A 2X2 matrix of two points in the second line defined the
%                same way as line1.
%
%OUTPUTS:   point The 2X1 intersection point of the two lines given as
%                 [x;y] components.
%
%The formula in terms of matrix determinants is taken from [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Line-Line Intersection." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/Line-LineIntersection.html
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

detL1=det(line1);
detL2=det(line2);

deltaL1=-diff(line1,1,2);
deltaL2=-diff(line2,1,2);

denom=det([deltaL1,deltaL2]);

x=det([[detL1;detL2],[deltaL1(1);deltaL2(1)]])/denom;
y=det([[detL1;detL2],[deltaL1(2);deltaL2(2)]])/denom;

point=[x;y];

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
