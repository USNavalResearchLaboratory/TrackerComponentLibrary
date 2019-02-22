function muCentral=raw2CentralMoments(muRaw)
%%RAW2CENTRALMOMENTS Given the raw (noncentral) moments of a scalar
%                    probability distribution, compute the central
%                    moments of the distribution.  
%
%INPUTS: muRaw An nX1 or 1Xn vector of scalar raw (noncentral) moments.
%              muRaw(n) is the nth raw moment. Thus, muRaw(1) is the mean
%              of the distribution. n>=1.
%
%OUTPUTS: muCentral An nX1 vector of the scalar central moments.
%                   muCentral(n) is the nth order central moment. Thus,
%                   muCentral(1)=0.
%
%The conversion between central and non-central moments is taken from
%Chapter 5.4 of [1]. Note that this function is implemented to work with
%exact raw moments. If sample raw moments are given, then muCentral will be
%biased. If there are not too many samples, then the function
%unbiasedMomentCumulant can be used to obtain unbiased sample central
%moments directly.
%
%The nth noncentral moment of the distribution is the expected value
%E{x^n}. The nth central moment of the distribution is E{(x-E\{x\})^n}.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(muRaw);

muCentral=zeros(n,1);
%Precompute powers to avoid having to recompute them in the loop.
mu1Pows=muRaw(1).^(0:n);

muCentral(1)=0;
for i=2:n
    %The zeroth-order term
    binomCoeff=1;
    signVal=-1*mod(i,2)+mod(i+1,2);
    muCentral(i)=signVal*binomCoeff*mu1Pows(i+1);

    for k=1:i
        binomCoeff=((i-k+1)/k)*binomCoeff;
        signVal=-signVal;
        
        muCentral(i)=muCentral(i)+signVal*binomCoeff*muRaw(k)*mu1Pows(i-k+1);
    end
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
