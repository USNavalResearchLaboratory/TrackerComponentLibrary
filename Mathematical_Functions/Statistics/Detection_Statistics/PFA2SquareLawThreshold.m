function thresh=PFA2SquareLawThreshold(PFA,N,ampDef)
%%PFA2SQUARELAWTHRESHOLD Given a desired probability of false alarm (PFA),
%        obtain the normalized (assuming average noise variance is 1 or 2
%        depending on the normalization) detection threshold for a square
%        law detector that would result in the PFA. This assumes corruption
%        with complex Gaussian noise whose normalization can vary as
%        specified by the ampDef input. The square law detector for N
%        samples is
%        y=sum_{i=1}^N abs(r_i)^2
%        where y is compared to a threshold. This is the inverse of
%        threshold2SquareLawPFA.
%
%INPUTS: PFA The desired probabily of false alarm. A matrix can be passed
%          if multiple thresholds are desired at once, 0<=PFA<=1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). If this parameter is
%          omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see the comments at the top of
%          NonFlucSqLawD for detailed definitions). Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: thresh The normalized detection threshold(s) needed to achieve
%                the desired false alarm rate(s).
%
%Equation 10.4.20 in Chapter 10.4 of [1] relates the false alarm rate of a
%square law detector to the threshold. This function inverts that
%expression using an inverse incomplete gamma function. However, if
%ampDef=0, it also transforms the threshold to account for Equation 10.4-4,
%which defines the threshold in terms of a scaled square law detector.
%
%EXAMPLE 1:
%In this example, we generate samples of the noise according to ampDef=0
%and then get a PFA based on a provied threshold. We then give that PFA to
%this function to recover the originally provided threshold to show that it
%is correct.
% numSamples=1e6;
% ampDef=0;
% N=3;
% thresh=6;
% noise=(randn(N,numSamples)+1j*randn(N,numSamples));
% y=sum(abs(noise).^2,1);
% PFA=mean(y>thresh)
% threshBack=PFA2SquareLawThreshold(PFA,N,ampDef)
%
%EXAMPLE 2:
%This is the mostly the same as Example 1, except we are using ampDef=1 and 
%thus must normalize the noise samples differently.
% numSamples=1e6;
% ampDef=1;
% N=10;
% thresh=15;
% noise=(1/sqrt(2))*(randn(N,numSamples)+1j*randn(N,numSamples));
% y=sum(abs(noise).^2,1);
% PFA=mean(y>thresh)
% threshBack=PFA2SquareLawThreshold(PFA,N,ampDef)
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(ampDef))
    ampDef=1;
end

if(nargin<2||isempty(N))
	N=1; 
end

x=1-PFA;

%Inverting Equation 10.4.20.
thresh=gammaincinv(x,N);

if(ampDef==0)
    %Account for Equation 10.4-4 in [1] using a halved version of the
    %square law.
    thresh=thresh*2;
end
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
