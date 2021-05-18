function area=hypersphereArea(r,numDim)
%%HYPERSPHEREAREA Determine the surface area of a hypersphere of a given
%                 radius in a given number of dimensions.
%
%INPUTS:     r The radius of the hypersphere. r>0.
%       numDim The number of dimensions of the hypersphere. numDim>=1. One
%              dimensions is a line, two a circle, three a sphere, etc.
%
%OUTPUTS: area The area of the hypersphere.
%
%The surface area of a hypersphere is given in terms of the gamma function
%in [1].
%
%REFERENCES:
%[1] S. Li, "Concise formulas for the area and volume of a hyperspherical
%    cap," Asian Journal of Mathematics and Statistics, vol. 4, no. 1, pp.
%    66-70, 2011.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

area=2*pi^(numDim/2)/gamma(numDim/2)*r^(numDim-1);

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
