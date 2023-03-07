function [xi,w]=thirdOrderTetrahedronCubPoints()
%%THIRDORDERTETRAHEDRONCUBPOINTS Obtain third-order cubature points
%   for integration over a tetrahedron in 3D. The points and weights are
%   for the tetrahedron with vertices (1,0,0), (0,1,0), (0,0,1), and
%   (0,0,0), but can be transformed to any tetrahedron using
%   transformSimplexTetrahedronPts.
%
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard
%            tetrahedron.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard tetrahedron (1/6).
%
%This function implements the points given in [1] (8 points).
%
%EXAMPLE:
%Given the vertices of the simplex, we compare a third-order moment
%computed using these cubature points to one computed using
%monomialIntSimplex. The results are the same within typical finite
%precision limits.
% [xi,w]=thirdOrderTetrahedronCubPoints();
% alpha=[1;1;1];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntSimplex(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[-0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   -0.9689798150982901207378151874892027072,   0.18162379004944980942342872025562069427;
   -0.34367339496723662642072827083693243093,   -0.9689798150982901207378151874892027072,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427;
    -0.9689798150982901207378151874892027072,  -0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427;
   -0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,  -0.34367339496723662642072827083693243093,   0.18162379004944980942342872025562069427;
   -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.35171650060942837529461966476512015033,   0.15170954328388352390990461307771263906;
   -0.78390550020314279176487322158837338344,   0.35171650060942837529461966476512015033,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906;
    0.35171650060942837529461966476512015033,  -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906;
   -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,  -0.78390550020314279176487322158837338344,   0.15170954328388352390990461307771263906];
w=M(:,4);
xi=M(:,1:3)';
%Transform the points to the standard tetrahedron.
v1=[-1, 1,-1,-1;
    -1,-1,-1, 1;
    -1,-1, 1, -1];
v2=[1,0,0,0;
    0,1,0,0;
    0,0,1,0];
[A,d]=affineTransBetweenTetrahedra(v1,v2);
xi=bsxfun(@plus,A*xi,d);
w=w/8;
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
