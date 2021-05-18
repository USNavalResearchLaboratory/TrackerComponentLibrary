function KVals=EpanechnikovKernel(x)
%%EPANECHNIKOVKERNEL Evaluate the multivariate Epanechnikov kernel at one
%                    or more points. The points should have already been
%                    scaled by the inverse badnwidth. This Kernel has
%                    certain optimality properties in kernel density
%                    estimation algorithms.
%
%INPUTS: x An xDimXNs set of points.
%
%OUTPUTS: KVals A 1XNs set of values of the Epanechnikov kernel evluated at
%               the given points.
%
%The multivariate Epanechnikov kernel is given in Equation 4.4 of Chapter
%4.2.1 of [1]. The Kernel is equivalent to the PDF of the multivariate
%Pearson type II distribution with a=2. Thus, this function is just a
%wrapper for PearsonIID.PDF. 
%
%REFERENCES:
%[1] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

KVals=PearsonIID.PDF(x,2);
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
