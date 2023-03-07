function [xi,w]=seventhOrderTriangleCubPoints()
%%SEVENTHORDERTRIANGLECUBPOINTS Obtain seventh-order cubature points for
%   integration over a triangle in 2D. The points and weights are for the
%   triangle with vertices (1,0), (0,1), (0,0), but can be transformed to
%   any triangle using transformSimplexTriPoints.
%
%INPUTS: None
%
%OUTPUTS: xi A 2XnumCubPoints set of points for the standard triangle.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the triangle (1/2).
%
%This function implements the points given in [1] (15 points).
%
%EXAMPLE:
%Given the vertices of the simplex, we compare a seventh-order moment
%computed using these cubature points to one computed using
%monomialIntSimplex. The results are the same within typical finite
%precision limits.
% [xi,w]=seventhOrderTriangleCubPoints();
% alpha=[4;3];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntSimplex(alpha)
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[-0.93253870289082430257005654739836752385,     0.8650774057816486051401130947967350477,   0.033090100221584262071955896945834891238;
     0.8650774057816486051401130947967350477,   -0.93253870289082430257005654739836752385,   0.033090100221584262071955896945834891238;
   -0.93253870289082430257005654739836752385,   -0.93253870289082430257005654739836752385,   0.033090100221584262071955896945834891238;
   -0.51684523480919288209962646032443600082,   0.033690469618385764199252920648872001638,    0.25588834246031114556580247036929262133;
   0.033690469618385764199252920648872001638,   -0.51684523480919288209962646032443600082,    0.25588834246031114556580247036929262133;
   -0.51684523480919288209962646032443600082,   -0.51684523480919288209962646032443600082,    0.25588834246031114556580247036929262133;
  -0.051380614990563531580838528101628435036,   -0.89723877001887293683832294379674312993,    0.15417329237197213566964304166748277293;
   -0.89723877001887293683832294379674312993,  -0.051380614990563531580838528101628435036,    0.15417329237197213566964304166748277293;
  -0.051380614990563531580838528101628435036,  -0.051380614990563531580838528101628435036,    0.15417329237197213566964304166748277293;
   -0.90592671069480953331718004928623004897,    0.50856008110010635471247864925623992738,    0.11175746580639956167963262884202819059;
    0.50856008110010635471247864925623992738,   -0.90592671069480953331718004928623004897,    0.11175746580639956167963262884202819059;
   -0.60263337040529682139529859997000987842,    0.50856008110010635471247864925623992738,    0.11175746580639956167963262884202819059;
    0.50856008110010635471247864925623992738,   -0.60263337040529682139529859997000987842,    0.11175746580639956167963262884202819059;
   -0.60263337040529682139529859997000987842,   -0.90592671069480953331718004928623004897,    0.11175746580639956167963262884202819059;
   -0.90592671069480953331718004928623004897,   -0.60263337040529682139529859997000987842,    0.11175746580639956167963262884202819059];

 w=M(:,3);
xi=M(:,1:2)';
%Transform the points to the standard triangle.
v1=[-1,-1, 1;
    -1, 1,-1];
v2=[1,0,0;
    0,1,0];
[A,d]=affineTransBetweenTriangles(v1,v2);
xi=bsxfun(@plus,A*xi,d);
w=w/4;
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
