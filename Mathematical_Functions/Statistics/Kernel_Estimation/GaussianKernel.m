function KVals=GaussianKernel(x)
%%GAUSSIANKERNEL Evaluate the multivariate Gaussian kernel at one or more
%                points. The points should have already been scaled by the
%                inverse badnwidth. This is just the zero-mean multivariate
%                Gaussian PDF with the identity matrix as its covariance
%                matrix.
%
%INPUTS: x An xDimXNs set of points.
%
%OUTPUTS: KVals A 1XNs set of values of the Gaussian kernel evluated at
%               the given points.
%
%Gaussian kernels play a role in the regularization step of particle
%filters, as mentioned in [1].
%
%This function just directly implements the formula, since a call to the
%PDF function in  GaussianD function would be slower (since it allowd a 
%
%REFERENCES:
%[1] C. Musso, N. Oudjane, and F. LeGland, "Improving regularised particle
%    filters," in Sequential Monte Carlo Methods in Practice, A. Doucet, J.
%    F. G. de Freitas, and N. J. Gordon, Eds., Ch. 12, New York:
%    Springer-Verlag, 2001.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);
KVals=(2*pi)^(-numDim/2)*exp(-1/2*sum(x.*x,1));

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
