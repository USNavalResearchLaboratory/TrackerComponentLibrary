function PFA=threshold2SquareLawPFA(thresh,N,ampDef)
%%THRESHOLD2SQUARELAWPFA Given a normalized detection threshold obtain the
%        probability of false alarm for a square law detector. This assumes
%        corruption with complex Gaussian noise whose normalization can
%        vary as specified by the ampDef input. The square law detector for
%        N samples is
%        y=sum_{i=1}^N abs(r_i)^2
%        where y is compared to a threshold. This is the inverse of
%        PFA2SquareLawThreshold.
%
%INPUTS: thresh The scalar normalized detection threshold to use. This is
%          the threshold to use if the noise variance is 1. A matrix can be
%          passd if multiple PFAs are desired in parallel.
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
%OUTPUTS: PFA The probability of false alarm, or a matrix of probabilities
%             if a matrix is passed for thresh, 0<=PFA<=1.
%
%Equation 10.4.20 in Chapter 10.4 of [1] relates the false alarm rate of a
%square law detector to the threshold.  However, if
%ampDef=0, the threshold has to be first transformed to reflect the
%different normalization used.
%
%EXAMPLE 1:
%Here, we verify that the false alarm rate predicted by this function is
%the same as that obtained through Monte Carlo runs using ampDef=0. We
%choose a threshold for a high PFA so that we do not need importance
%sampling to get a good result. The computed and Monto Carlo PFAs will be
%close.
% numSamples=1e6;
% N=4;
% thresh=6;
% ampDef=0;
% PFAComputed=threshold2SquareLawPFA(thresh,N,ampDef)
% noise=(randn(N,numSamples)+1j*randn(N,numSamples));
% y=sum(abs(noise).^2,1);
% PFAMonteCarlo=mean(y>=thresh)
%
%EXAMPLE 2:
%This is the same as example 1, except ampDef=1 and the generation of the
%noise matches the new definition.
% numSamples=1e6;
% N=4;
% thresh=6;
% ampDef=1;
% PFAComputed=threshold2SquareLawPFA(thresh,N,ampDef)
% noise=(1/sqrt(2))*(randn(N,numSamples)+1j*randn(N,numSamples));
% y=sum(abs(noise).^2,1);
% PFAMonteCarlo=mean(y>=thresh)
%
%EXMAPLE 3:
%THis computed a PFA from a threshold and then uses PFA2SquareLawThreshold
%to reconstruct the threshold, The relative error of the inverted value is
%in the order of finite precision errors. The relative error can increase
%as the threshold becomes very large.
% N=10;
% thresh=30;
% ampDef=1;
% PFA=threshold2SquareLawPFA(thresh,N,ampDef);
% threshBack=PFA2SquareLawThreshold(PFA,N,ampDef);
% RelErr=(thresh-threshBack)/thresh
%
%REFERENCES:
%[1]J. V. Di Franco and W. L. Rubin, Radar Detection. Prentice Hall, Inc.,
%   Englewood: 1968.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(ampDef))
    ampDef=1;
end

if(nargin<2||isempty(N))
	N=1; 
end

if(ampDef==0)
    thresh=thresh/2;
end

%This is Equation 10.4.20 in Chapter 10.4 of [1].
PFA=1-PearsonsGammaInc(thresh./sqrt(N),N-1);

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
