function [xi,w]=sixthOrderTriangleCubPoints()
%%SIXTHORDERTRIANGLECUBPOINTS Obtain sixth-order cubature points for
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
%This function implements the points given in [1] (12 points).
%
%EXAMPLE:
%Given the vertices of the simplex, we compare a sixth-order moment
%computed using these cubature points to one computed using
%monomialIntSimplex. The results are the same within typical finite
%precision limits.
% [xi,w]=sixthOrderTriangleCubPoints();
% alpha=[4;2];
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

M=[-0.87382197101699554331933679425836168532,    0.74764394203399108663867358851672337064,    0.10168981274041363384187361821373796809;
    0.74764394203399108663867358851672337064,   -0.87382197101699554331933679425836168532,    0.10168981274041363384187361821373796809;
   -0.87382197101699554331933679425836168532,   -0.87382197101699554331933679425836168532,    0.10168981274041363384187361821373796809;
   -0.50142650965817915741672289378596184782,  0.0028530193163583148334457875719236956481,    0.23357255145275873205057922277115888265;
  0.0028530193163583148334457875719236956481,   -0.50142650965817915741672289378596184782,    0.23357255145275873205057922277115888265;
   -0.50142650965817915741672289378596184782,   -0.50142650965817915741672289378596184782,    0.23357255145275873205057922277115888265;
   -0.89370990031036610529350065673720370601,    0.27300499824279729446028518882409939961,    0.16570215123674715038710691284088490796;
    0.27300499824279729446028518882409939961,   -0.89370990031036610529350065673720370601,    0.16570215123674715038710691284088490796;
   -0.37929509793243118916678453208689569359,    0.27300499824279729446028518882409939961,    0.16570215123674715038710691284088490796;
    0.27300499824279729446028518882409939961,   -0.37929509793243118916678453208689569359,    0.16570215123674715038710691284088490796;
   -0.37929509793243118916678453208689569359,   -0.89370990031036610529350065673720370601,    0.16570215123674715038710691284088490796;
   -0.89370990031036610529350065673720370601,   -0.37929509793243118916678453208689569359,    0.16570215123674715038710691284088490796];
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
