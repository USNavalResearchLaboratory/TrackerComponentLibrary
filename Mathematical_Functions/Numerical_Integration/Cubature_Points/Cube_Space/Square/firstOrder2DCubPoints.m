function [xi,w]=firstOrder2DCubPoints()
%%FIRSTORDER2DCUBPOINTS Generate first-order cubature points
%               for integration over a 2-dimensional cube with bounds of 
%               (-1,-1), (-1,1), (1,1), and (1,-1).
%
%INPUTS: None
%
%OUTPUTS: xi This is a 2XnumCubPoints set of points for the standard
%            square.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard square (4).
%
%This function implements the points given in [1] (1 point). It is just the
%origin with a weight of 4.
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

w=4;
xi=[0;0;0];

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
