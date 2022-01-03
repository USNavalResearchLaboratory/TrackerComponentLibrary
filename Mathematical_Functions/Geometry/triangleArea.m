function [A,dAdv1]=triangleArea(v1,v2,v3,onlyPositive)
%%TRIANGLEAREA Find the signed area of a triangle in 2D given 3 points. The
%          sign is positive if the vertices are given counterclockwise and
%          it is negative if they are given clockwise.
%
%INPUTS: v1,v2,v3 A 2X1 vectors in order of the edges of the triangle.
%    onlyPositive An optional parameter. If true, then only non-negative
%                 values of A will be returned (an unsigned area).
%                 Otherwise, A is signed. The default if omitted or an
%                 empty matrix is passed is false.
%
%OUTPUTS: A The signed area of the triangle (or positive value if
%           onlyPositive is true).
%     dAdv1 The 1X2 derivative of the signed area with respect to the
%           elements of v1 in the order dAdv1=[dA/dv1(1),dA/dv1(2)].
%
%The signed area of a triangle can be expressed using a determinant as in
%[1]. The sum is implemented here using accurateSum.
%
%EXAMPLE 1:
%Here, the signed area is -0.5.
% v1=[0;
%     0];
% v2=[0;
%     1];
% v3=[1;
%     1];
% A=triangleArea(v1,v2,v3)
%
%EXAMPLE 2:
%In this example, we verify that the derivative is consistent with
%numerical differentiation. The absolute error is consistent with what one
%might expect due to finite precision limitations.
% v1=[0;
%     0];
% v2=[0;
%     1];
% v3=[1;
%     1];
% [~,dAdx]=triangleArea(v1,v2,v3);
% f=@(x)triangleArea(x,v2,v3);
% dadxNumDiff=numDiff(v1,f,1);
% AbsErr=dadxNumDiff-dAdx
%
%REFERENCES:
%[1] Weisstein, Eric W. "Triangle Area." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/TriangleArea.html
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(onlyPositive))
    onlyPositive=false;
end

x1=v1(1);
x2=v2(1);
x3=v3(1);
y1=v1(2);
y2=v2(2);
y3=v3(2);

A=accurateSum([x1*y2;-x2*y1;x3*y1;-x3*y2;-x1*y3;x2*y3],0,true)/2;

if(onlyPositive)
    s=sign(A);
    A=s*A;
else
    s=1;
end

if(nargout>1)
    dAdv1=s*[y2-y3,x3-x2]/2;
end
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
