function vals=sampleEpanechnikovKernel(N,d)
%%SAMPLEEPANECHNIKOVKERNEL Draw samples of the Epanechnikov kernel. This
%           can be useful when performing regularization in particle
%           filters, as described in Chapter 3.5.3 of [1].
%
%INPUTS: N The number of random samples to draw.
%        d The dimensionality of the samples.
%
%OUTPUTS: vals The dXN set of N samples.
%
%The Epanechnikov kernel is equivalent to the multivariate Pearson type II
%distribution with a=2. Thus, this function is just a wrapper for
%PearsonIID.rand.
%
%REFERENCES:
%[1] B. Ristic, S. Arulampalam, and N. Gordon, Beyond the Kalman Filter:
%    Particle Filters for Tracking Applications. Boston: Artech House,
%    2004.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

vals=PearsonIID.rand(N,d,2);

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
