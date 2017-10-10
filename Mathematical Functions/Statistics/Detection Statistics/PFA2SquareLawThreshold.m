function thresh=PFA2SquareLawThreshold(PFA,N)
%%PFA2SQUARELAWTHRESHOLD Given a desired probability of false alarm (PFA),
%        obtain the normalized (assuming average noise variance is 1)
%        detection threshold for a square law detector that would result in
%        the PFA. This assumes corruption with Gaussian noise with unit
%        variance.
%
%INPUTS: PFA The desired probabily of false alarm. A matrix can be passed
%            if multiple thresholds are desired at once, 0<=PFA<=1.
%          N The number of pulses that are to be incoherently added for
%            detection (in a square-law detector). If this parameter is
%            omitted or an empty matrix is passed, N=1 is used.
%
%OUTPUTS: thresh The normalized detection threshold(s) needed to achieve
%                the desired false alarm rate(s).
%
%The square law detector for N samples is
% y=sum_{i=1}^N r_i^2
%where r_i is the ith real amplitude sampled. The model comes from taking a
%sample of a complex amplitude coming from a filter. The squared real
%amplitude is r^2_i=y_{I,i}^2+y_{Q,i}^2 where y_{I,i} and y_{Q,i} are the
%in-phase and quadrature components of the filter output. In the absence of
%a signal, the value y=y_{I,i}+sqrt(-1)*y_{Q,i} in the absence of a signal
%is modeled as being distributed circularly complex Gaussian with zero mean
%and variance 2. (See, Equation 9.3-35a in [1]). The variance being 2
%simply reflects having normalized the variance on the I and Q components
%each to 1, which means that the threshold must be similarly normalized for
%use in this function.
%
%Equation 10.4.20 in Chapter 10.4 of [1] relates the false alarm rate of a
%square law detector to the threshold. This function inverts that
%expression using an inverse incomplete gamma function. However, it
%transforms the threshold to account for Equation 10.4-4, which defines the
%threshold in terms of a scaled square law detector.
%
%Under this model, one can generate a random false alarm sample from the
%complex Gaussian distribution as
% sample=sum(abs(ComplexGaussianD.rand(N1,0,2)).^2);
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(N))
	N=1; 
end

x=1-PFA;

%Inverting Equation 10.4.20.
thresh=gammaincinv(x,N);

%Account for Equation 10.4-4 in [1] using a halved version of the square
%law.
thresh=thresh*2;

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
