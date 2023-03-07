function [xi,w]=transformSimplexTriPts(xi,w,v)
%%TRANSFORMSIMPLEXTRIPTSTOTRI Given cubature points and weights for
%   integration over a the standard simplex triangle with vertices (1,0),
%   (0,1), and (0,0), transform the points and weights for integration over
%   the triangle with vertices given in v. 
%
%INPUTS: xi The 2XnumCubPoints set of cubature points for the standard 2D
%           simplex triangle.
%         w The numCubPointsX1 associated cubature weights for the standard
%           2D simplex triangle.
%         v The 2X3 set of vertices of the new triangle.
%
%OUTPUTS: xi The 2XnumCubPoints set of cubature points for the specified
%            triangle.
%          w The 2XnumCubPoints set of transformed weights.
%
%With cubature formulae as in [1], the weights must sum to the area of the
%triangle. Thus, the first step here is to scale the weights by 2*the
%triangle area (the 2 is because the area of the basic simplex triangle is
%1/6). An affine transformation, as described in Chapter 5 of [2], is then
%performed to move the cubature weights into the new triangle.
% 
%EXAMPLE:
%In this example, we find the integral of x^3*y^2 using fifth-order simplex
%points. The exact solution is 128/5 and one can see that the relative
%error between the solution found and the exact solution is within what one
%would expect due to finite precision limits.
% vertices=[0,2,4;
%           0,2,0];
% [xi,w]=fifthOrderSimplexCubPoints(2);
% [xi,w]=transformSimplexTriPts(xi,w,vertices);
% alpha=[3;2];
% cubSol=findMomentFromSamp(alpha,xi,w)
% exactSol=128/5
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] K. Anjyo and H. Ochiai, Mathematical Basics of Motion and Deformation
%    in Computer Graphics, 2nd ed. Morgan and Claypool Publishers, 2017.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Adjust w for the triangle area.
T=triangleArea(v(:,1),v(:,2),v(:,3),true);
w=w*(2*T);

%Transforms the points into the new triangle.
A=[v(:,1)-v(:,3),v(:,2)-v(:,3)];
xi=A*xi+v(:,3);

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
