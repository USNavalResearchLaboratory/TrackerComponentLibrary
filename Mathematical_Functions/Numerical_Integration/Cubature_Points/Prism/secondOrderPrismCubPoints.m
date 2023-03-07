function [xi,w]=secondOrderPrismCubPoints()
%%SECONDORDERPRISMCUBPOINTS Generate second-order cubature points for
%  integration over a standard prism in 3D. This is a triangle that has
%  been extruded upward. The vertices are (1,-1,-1), (-1,-1,-1), (-1,1,-1),
%  (1,-1,1), (-1,-1,1), and (-1,1,1).
% 
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard prism.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard prism (4).
%
%This function implements the points given in [1] (5 points).
%
%EXAMPLE:
%We compare a 2nd-order moment computed using these cubature points
%to one computed using monomialIntPrism. The results are the same within
%typical finite precision limits.
% [xi,w]=secondOrderPrismCubPoints();
% alpha=[2;0;0];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntPrism(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[-0.33333333333333333333333333333333333333,  -0.33333333333333333333333333333333333333,  -0.99999999999999985119303811076522004183,   0.66666666666666686507594918564641756458;
-0.33333333333333333333333333333333333333,  -0.33333333333333333333333333333333333333,   0.99999999999999985119303811076522004183,   0.66666666666666686507594918564641756458;
-0.74158162379719638007464124598498266439,   0.48316324759439276014928249196996532878,                                          0,   0.88888888888888875661603387623572162361;
 0.48316324759439276014928249196996532878,  -0.74158162379719638007464124598498266439,                                          0,   0.88888888888888875661603387623572162361;
-0.74158162379719638007464124598498266439,  -0.74158162379719638007464124598498266439,                                          0,   0.88888888888888875661603387623572162361];

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
