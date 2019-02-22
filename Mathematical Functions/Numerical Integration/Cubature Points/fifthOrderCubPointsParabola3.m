function [xi,w]=fifthOrderCubPointsParabola3(a,b)
%%FIFTHORDERCUBPOINTSPARABOLA3 Generate ffith order cubature points for
%           performing second order integration of a function over the
%           region between two parabolas with equation y=b-b*x^2/a^2 and
%           y=b*x^2/a^2-b. 
%
%INPUTS: a,b The two real parameters of the parabolas, whose equations are
%            given above.
%
%OUTPUTS: xi A 2XnumCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%THis function implements Formula Par^3 5-1 in [1], pg. 338, 13 points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xi=[[0;0],PMCombos([a/2;0]),PMCombos([a;0]),PMCombos([0;b/2]),PMCombos([0;b]),PMCombos([a/2;b/2])];

V=8*a*b/3;

B0=(344/6930)*V;
B1=(704/6930)*V;
B2=(165/6930)*V;
B3=(768/6930)*V;
B4=(248/6930)*V;
B5=B1;

w=[B0;B1;B1;B2;B2;B3;B3;B4;B4;B5*ones(4,1)];

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
