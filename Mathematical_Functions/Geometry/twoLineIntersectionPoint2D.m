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
%The formula in terms of matrix determinants is taken from [1]. It was
%observed that certain orderings of the terms can lead to poor finite
%precision performance.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Line-Line Intersection." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/Line-LineIntersection.html
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

x1=line1(1,1);
y1=line1(2,1);
x2=line1(1,2);
y2=line1(2,2);
x3=line2(1,1);
y3=line2(2,1);
x4=line2(1,2);
y4=line2(2,2);

num=(x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4); 
den=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);

x=num/den;

num=(x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4); 
den=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4); 
y=num/den;
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
