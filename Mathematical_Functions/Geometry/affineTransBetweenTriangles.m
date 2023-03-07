function [A,d]=affineTransBetweenTriangles(v1,v2)
%%AFFINETRANSBETWEENTRIANGLES Given the vertices of one triangle in 2D and
%       the vertices of a second triangle in 2D, obtain a matrix M and a
%       vector d such that p2=M*p1+d transforms a point p1 in the first
%       triangle into a point p2 in the second triangle.
%
%INPUTS: v1 The 2X3 set of vertices in the first triangle.
%        v2 The 2X3 set of vertices in the second triangle (into which we
%           wish to transform the first triangle).
%
%OUTPUTS: A A 2X2 scale/rotation matrix.
%         d A 2X1 origin shift vector.
%
%Equations for the transformation are given in Chapter 5 of [1].
%
%EXAMPLE:
%This just find the transformation between two sets of vertices and
%transforms the first set of vertices into the second set. The error is on
%the order of typical finite precision limitations.
% v1=[ 5, -22,   3;
%     18,   8, -13];
% v2=[-4, 35, -13;
%      3, 27,  30];
% [A,d]=affineTransBetweenTriangles(v1,v2);
% v2Trans=bsxfun(@plus,A*v1,d);
% RelErr=max(abs((v2Trans-v2)./v2))
%
%REFERENCES:
%[1] K. Anjyo and H. Ochiai, Mathematical Basics of Motion and Deformation
%    in Computer Graphics, 2nd ed. Morgan and Claypool Publishers, 2017.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A1=[v1;
   1,1,1];
A2=[v2;
   1,1,1];

M=A2/A1;

A=M(1:2,1:2);
d=M(1:2,3);

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
