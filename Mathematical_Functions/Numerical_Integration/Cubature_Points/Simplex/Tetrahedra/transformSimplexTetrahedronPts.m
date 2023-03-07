function [xi,w]=transformSimplexTetrahedronPts(xi,w,v)
%%TRANSFORMSIMPLEXTETRAHEDRONPTS Given cubature points and weights for
%   integration over a the standard simplex tetrahedron with vertices
%   (1,0,0), (0,1,0), (0,0,1) and (0,0,0), transform the points and weights
%   for integration over the tetrahedron with vertices given in v. 
%
%INPUTS: xi The 3XnumCubPoints set of cubature points for the standard 2D
%           simplex tetrahedron.
%         w The numCubPointsX1 associated cubature weights for the standard
%           3D simplex tetrahedron.
%         v The 3X4 set of vertices of the new tetrahedron.
%
%OUTPUTS: xi The 3XnumCubPoints set of cubature points for the specified
%            tetrahedron.
%          w The 3XnumCubPoints set of transformed weights.
%
%With cubature formulae as in [1], the weights must sum to the area of the
%tetrahedron. Thus, the first step here is to scale the weights by 6*the
%tetrahedron area (the 6 is because the area of the basic simplex
%tetrahedron is 1/6). A 3D affine transformation, analogous to that used
%for triangles in 2D Chapter 5 of [2], is then performed to move the
%cubature weights into the new tetrahedron.
%
%EXAMPLE:
%As an example, we take the tetrahedron with vertices (2,0,0), (0,3,0),
%(0,0,1) and (0,0,0). The triple integral over this region is of the form:
%int_0^1 int_0^{2(1-z)} int_0^{3(1-x/2-z)} dy dx dz
%We consider an eight order moment where the argument of the integral is
%x^4*y^2*z^2 The exact analytic solution is 4/1925. We compare that to the
%solution using transformed eigth order cubature points and we see that
%the relative error is within typical finite precision limits, so this
%transformation works.
% vertices=[2,0,0,0;
%           0,3,0,0;
%           0,0,1,0];
% [xi,w]=eighthOrderTetrahedronPoints();
% [xi,w]=transformSimplexTetrahedronPts(xi,w,vertices);
% alpha=[4;2;2];
% cubSol=findMomentFromSamp(alpha,xi,w);
% exactSol=4/1925;
% RelErr=(cubSol-exactSol)/exactSol
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] K. Anjyo and H. Ochiai, Mathematical Basics of Motion and Deformation
%    in Computer Graphics, 2nd ed. Morgan and Claypool Publishers, 2017.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Adjust w for the tetrahedron volume
T=tetrahedronVolume(v,true);
w=w*(6*T);

%Transforms the points into the new triangle.
A=[v(:,1)-v(:,4),v(:,2)-v(:,4),v(:,3)-v(:,4)];
xi=bsxfun(@plus,A*xi,v(:,4));

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
