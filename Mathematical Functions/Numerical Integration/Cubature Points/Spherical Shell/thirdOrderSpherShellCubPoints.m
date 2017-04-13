function [xi,w]=thirdOrderSpherShellCubPoints(numDim,rMin)
%%THIRDORDERSPHERSHELLCUBPOINTS Obtain third order cubature points and
%               weight for integrating spherical region surrounding the
%               origin from radius rMin to radius 1, where rMin>0.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%          rMin The minimum radius of the spherical shell. 0<rMin<1.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This implements formula S_n^{shell} 3-1 in [1], pg. 293, 2*numDim points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

rho=rMin;
n=numDim;
V=(2*pi^(n/2)/(n*gamma(n/2)))*(1-rho^n);

r=sqrt(n/(n+2)*((1-rho^(n+2))/(1-rho^n)));
B=V/(2*n);

xi=fullSymPerms([r;zeros(n-1,1)]);
w=B*ones(2*n,1);

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
