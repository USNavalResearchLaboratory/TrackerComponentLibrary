function [xi,w]=secondOrderPyramidCubPoints()
%%SECONDORDERPYRAMIDCUBPOINTS Generate second-order cubature points for
%  integration over a 3-dimensional pyramid with a square base with the
%  peak vertex at (0,0,1) and the base vertices at (1,-1,-1), (-1,-1,-1),
%  (-1,1,-1), and (1,1,-1).
% 
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard
%            square pyramid.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard square pyramid (8/3).
%
%This function implements the points given in [1] (5 points).
%
%EXAMPLE:
%We compare a 2nd-order moment computed using these cubature points
%to one computed using monomialIntPyramid. The results are the same within
%typical finite precision limits.
% [xi,w]=secondOrderPyramidCubPoints();
% alpha=[0;2;0];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntPyramid(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[                                     0,                                          0,   0.21658207711955775339238838942231815011,   0.60287280353093925911579186632475611728;
 0.71892105581179616210276971993495767914,                                          0,  -0.70932703285428855378369530000365161136,   0.51594846578393185188771870008547763735;
                                        0,   0.71892105581179616210276971993495767914,  -0.70932703285428855378369530000365161136,   0.51594846578393185188771870008547763735;
-0.71892105581179616210276971993495767914,                                          0,  -0.70932703285428855378369530000365161136,   0.51594846578393185188771870008547763735;
                                        0,  -0.71892105581179616210276971993495767914,  -0.70932703285428855378369530000365161136,   0.51594846578393185188771870008547763735];

w=M(:,4);
xi=M(:,1:3)';

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
