function val=ItoIntegral2ndMoment(alphaIdx,betaIdx,deltaT)
%%ITOINTEGRAL2NDMOMENT Find the expected value of the product of two Ito
%                     integrals with multiple indices, specifically,
%                     val=E[I_{alphaIdx}*I_{betaIdx}]
%                     where I_{alphaIdx} is multiple Ito integrals of an
%                     argument of 1 from 0 to deltaT. The indices indicate
%                     which of multiple indices in a vector Wiener process
%                     is the measure for a given integral. An index of 0
%                     indicates that that integral is deterministic (dt not
%                     dW(i) where i is some index) This function only finds
%                     solutions where the nonzero elements of alphaIdx and
%                     betaIdx are equal and in the same order (though there
%                     can be a varying order of 0 terms between elements),
%                     which means that this function finds second moments,
%                     but not cross terms.
%
%INPUTS: alphaIdx The mX1 or 1Xm vector of indices of integration for the
%                 first Ito integral. The indices are integers >=0. 0
%                 indicates that the integral is taken over 
%         betaIdx The indices of the second Ito integral. This has the same
%                 format as alphaIdx. This must have the same number of
%                 nonzero components as alphaIdx and they must be in the
%                 same order as in alphaIdx.
%          deltaT The timespan of integration. deltaT>=0.
%
%OUTPUTS: val The real, scalar value of the desired second noncentral
%             moment.
%
%This function implements Eqautions 7.6 and 7.7 in Chapter 5.7 of [1].
%
%EXAMPLE:
%Here, we compare the value obtained by this function to an exact solution
%for a particular Itô Integral. In Chapter 10.4, an explicit expression for
%I_{j,j,j} is given. It is 
%I_{j,j,j}=(1/2)*((1/3)*(sqrt(deltaT)*xi)^2-deltaT)*sqrt(deltaT)*xi
%where xi is a normal random variable with zero mean and unit variance.
%Squaring that and taking the expected value, one gets
%exactSol=deltaT^3*(m2/4-m4/6+m6/36)
%where m2, m4, and m6 and the second, fourth and sixth noncentral moment of
%the normal-0-1 distribution. These are m2=1, m4=3, m6=15. Thus, the
%solution is exactSol=deltaT^3/6.
% deltaT=12;
% exactSol=deltaT^3/6
% val=ItoIntegral2ndMoment([1,1,1],[1,1,1],deltaT)
%One will see that the value returned by this function agrees with the
%exact solution, which is 288.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[alphaPlus,kAlpha]=getKAlphaPlus(alphaIdx);
[betaPlus,kBeta]=getKAlphaPlus(betaIdx);

lAlpha=length(alphaPlus);
lBeta=length(betaPlus);

if(lAlpha~=lBeta)
    error('This function only works when alphaIdx and betaIdx have the same number of nonzero components.')
end 

if(any(alphaPlus~=betaPlus)||deltaT==0)
    val=0;%Equation 7.6
    return;
end

%Equation 7.5.
wab=lAlpha+sum(kAlpha+kBeta);

if(wab<=170)
    CProd=1;
    for i=1:(lAlpha+1)
        CProd=CProd*binomial(kAlpha(i)+kBeta(i),kAlpha(i));
    end

    val=deltaT^wab/factorial(wab)*CProd;
else
    %Use logarithms of values in intermediate steps when overflows might
    %occur.
    logCProd=0;
    for i=1:(lAlpha+1)
        n=kAlpha(i)+kBeta(i);
        k=kAlpha(i);
        logCProd=logCProd+gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1);
    end

    val=exp(wab*log(deltaT)-gammaln(wab+1)+logCProd);
end
end

function [alphaPlus,k]=getKAlphaPlus(alpha)
%%GETKALPHAPLUS Thus function can be used to find alpha+ and k(alpha) as
%               well as beta+ and k(beta) as described in Section 5.7 of
%               [1].
%
%alpha+ is just alpha without the zero elements. The first element in
%k(alpha) is the number of zero element prior to the first nonzero element.
%k(i) is then the number of zero elements between the ith and the i+1th
%nonzero component of alpha or the end of alpha.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=length(alpha);
k=zeros(n+1,1);
alphaPlus=zeros(n,1);
numPlus=0;

idx=1;
while(idx<=n&&alpha(idx)==0)
    idx=idx+1;
end
k(1)=idx-1;%The k0 term.

if(idx<=n)
    numPlus=numPlus+1;
    alphaPlus(numPlus)=alpha(idx);
else%It is all zeros.
    alphaPlus=[];
    k=k(1);
    return;
end

gapLength=0;
idx=idx+1;
while(idx<=n)
    if(alpha(idx)==0)
        gapLength=gapLength+1;
    else
        numPlus=numPlus+1;
        alphaPlus(numPlus)=alpha(idx);
        k(numPlus)=gapLength;
        gapLength=0;
    end
    idx=idx+1;
end

k(numPlus+1)=gapLength;
%Size to fit.
alphaPlus=alphaPlus(1:numPlus);
k=k(1:(numPlus+1));

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
