function muRaw=central2RawMoments(muCentral,mu1)
%%CENTRAL2RAWMOMENTS Given the central moments of a scalar probability
%                    distribution, compute the central moments of the
%                    distribution.
%
%INPUTS: muCentral An nX1 vector of the scalar central moments.
%                  muCentral(n) is the nth order central moment. Thus,
%                  muCentral(1)=0. n>=1.
%              mu1 The scalar mean of the distribution.
%
%OUTPUTS: muRaw An nX1 or 1Xn vector of scalar raw (noncentral) moments.
%               muRaw(n) is the nth raw moment. Thus, muRaw(1) is the mean
%               of the distribution.
%
%The conversion between central and non-central moments it aken from
%Chapter 5.4 of [1].
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

n=length(muCentral);

muRaw=zeros(n,1);

muRaw(1)=mu1;
%Precompute powers to avoid having to recompute them in the loop.
mu1Pows=mu1.^(0:n);
for i=2:n
    %The zeroth-order term
    binomCoeff=1;
    muRaw(i)=binomCoeff*mu1Pows(i+1);
    
    for k=1:i
        binomCoeff=((i-k+1)/k)*binomCoeff;
        
        muRaw(i)=muRaw(i)+binomCoeff*muCentral(k)*mu1Pows(i-k+1);
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
