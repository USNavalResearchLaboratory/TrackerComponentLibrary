function [coeffs,polesRep]=partialFracKnownPoleDenom(poles,k,numerPoly)
%%PARTIALFRACKNOWNPOLEDENOM Given the numerator and the factored
%  denominator of a fraction, numerator/denominator with
%  numerator=numerPoly(1)*x^nn+numerPoly(2)*x^(nn-1)....+numerPoly(nn+1)*x+numerPoly(nn+1)
%  and
%  denominator=(x-poles(1))^k(1)*(x-poles(2))^k(2)*...*(x-poles(nd))^k(nd)
%  Obtain the coefficients of the partial fraction expansion, the format of
%  which is discussed below. Unlike Matlab's residue function, having been
%  passed the explicit poles (the roots of the denominator polynomia)
%  allows the expansion to be computed to a much higher accuracy with many
%  repeated poles than if just a polynomial were passed. This function does
%  not allow the numerator to be of a higher polynomial order than the
%  denominator. Partial fractions are useful when performing inverse Z and
%  Laplace transforms.
% 
%INPUTS: poles An ndX1 or 1Xnd vector of the roots of the polynomial in the
%              denominator (the poles). Repeated roots should NOT be
%              repeated here. Rather, they should be labeled as being
%              repeated using the k input.
%            k An ndX1 or 1Xnd vector of the number of times each root is
%              repeated. Theses are integers >=1. If this is omitted or an
%              empty matrix is passed, it is assumed that none of the roots
%              are repeated.
%    numerPoly The nnX1 or 1Xnn coefficients of the polynomial in the
%              numerator given in non-factored form. This is in the same
%              format as is used in the polyval function. However, there
%              shouldn't be any leading zeros (i.e higher order terms that
%              are not used). The highest order coefficient comes first and
%              missing terms are given zero weight. If this is omitted or
%              an empty matrix is passed, then a value of 1 is used.
%
%OUTPUTS: coeffs A numAX1 vector of the coefficients of the vectors in the
%                partial fraction expansion. See the format below.
%       polesRep These are the numAX1 set of poles corresponding to each
%                coefficient. This will be the same as poles on the input
%                except the ith poles is repeated k(i) times.
%
%The partial fraction expansion lets one replace the above
%numerator/denominator fraction with a series of the form
%partialFrac=coeffs(1)/(x-poles(1))+coeffs(1)/(x-poles(1))^2+...coeffs(k(1))/(x-poles(1))^k(1)+
%coeffs(k(1)+1)/(x-poles(2))+...+coeffs(k(1)+k(2))/(x-poles(2))^k(2)+...
%coeffs(k(1)+k(2)+1)/(x-poles(3))+...+coeffs(k(1)+k(2)+k(3))/(x-poles(3))^k(3)...
%Thus, the ith pole leads to k(i) fractions with denominators of increasing
%order. We already know the denominators and that each item in poles gets
%repeated by the corresponding number of times in k, so we just need to
%solve or coeffs.
%
%The basic notion behind solving for coeffs is that if we write
%numerator/denominator=partialFrac
%we can multiply both sides by denominator. This gets rid of all of the
%denominators on the partialFrac side. We are then left with all of the
%coeffs terms times polynomials in x. The solution comes from realizing
%that for each order of x from x^n to x^1 to x^0 on the right, it might
%equal the same thing in the numerator on the left. Thus, we can write a
%linear system of equations to solve to perform the partial fraction
%expansion.
%
%EXAMPLE 1:
%We will do a partial fraction expansion of 
%1/((x-1)^2*(x+2)*(x-12)^3*(x-4))
%Using this function and compare it to an analytic solution. The relative
%error will imply more than 13 digits of precision, which is good.
% poles=[1,-2,12,4];
% k=[2,1,3,1];
% [coeffsRep,polesRep]=partialFracKnownPoleDenom(poles,k)
% coeffsExact=[1/43923;1/11979;... %(x-1) terms
%              1/148176;... %(x+2) term
%              34213/5142387712;-233/8348032;1/13552;... %(x-12) terms
%              -1/27648];%(x-4) term.
% RelErr=max(abs((coeffsRep-coeffsExact)./coeffsExact))
%
%EXAMPLE 2:
%This example has a lot of repeated zeros. That would be very difficult for
%a normal partial fraction expansion algorithm that would have to find the
%poles and might made bad mistakes. Thus, this compares the benefits of
%this function to using Matlab's residue function, which cannot take
%advantage of knowing a factored form of the denominator. The output of
%this function has a low absolute error though there is a loss of precision
%as shown by the relative error only implying a bit more than 3 significant
%digits for some terms. On the other hand, as of MatlabR2022b, the residue
%function produces complex poles and complex coefficients that are far
%removed from what they should be as shown by a high absolute error and a
%relative error on the order of 1e16.
% poles=[1,12];
% k=[12,8];
% [coeffsRep,polesRep]=partialFracKnownPoleDenom(poles,k);
% %Because we gave it the exact poles, polesRep holds the exact poles. Since
% %poles going into partialFracKnownPoleDenom is sorted, polesRep will be
% %sorted.
% %The exact solutions coefficients were computed and are:
% coeffsExact=[31824/61159090448414546291;%Beginning of (x-1) terms.
%              1768/505447028499293771;
%              1040/45949729863572161;
%              585/4177248169415651;
%              312/379749833583241;
%              156/34522712143931;
%              72/3138428376721;
%              30/285311670611;
%              120/285311670611;
%              36/25937424601;
%              8/2357947691;
%              1/214358881;%End of (x-1) terms.
%              -31824/61159090448414546291;%Beginning of (x-12) terms.
%              12376/5559917313492231481;
%              -4368/505447028499293771;
%              1365/45949729863572161;
%              -364/4177248169415651;
%              78/379749833583241;
%              -12/34522712143931;
%              1/3138428376721];
% 
% AbsErr=max(abs(coeffsRep-coeffsExact))
% RelErr=max(abs((coeffsRep-coeffsExact)./coeffsExact))
% 
% denomPoly=flip([429981696,-5446434816,31902253056,-114532392960,...
%                 281591686656,-501770322432, 668881736640,...
%                 -678854801760, 528900908545, -316781074380,...
%                 145284570882, -50612518012, 13252488687,...
%                 -2580145560, 369984732, -38674872, 2900463,...
%                 -151708, 5250, -108,1]);
% [coeffsMatlab,polesMatlab]=residue(1,denomPoly);
% %We use a STABLE sorting so that when we put the poles in the correct
% %order, it doesn't change the ordering within the repeated poles.
% [polesMatlab,idxList]=bubbleSort(polesMatlab);
% coeffsMatlab=coeffsMatlab(idxList);
% 
% RelErrorPolesMatlab=max(abs((polesMatlab-polesRep)./polesRep))
% AbsErrCoeffsMatlab=max(abs(coeffsMatlab-coeffsExact))
% RelErrorCoeffsMatlab=max(abs((coeffsMatlab-coeffsExact)./coeffsExact))
%
%August 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(numerPoly))
    numerPoly=1;
end

numPoles=length(poles);
if(nargin<2||isempty(k))
    k=ones(numPoles,1);
end

if(numPoles~=length(k))
    error('The number of poles needs to equal the length of k.')
end

if(any(k<=0)||any(~isreal(k)))
    error('The multiplicities k must be real, positive and not zero.')
end

numA=sum(k);%The number of coefficients to find.
numNumerTerms=length(numerPoly);
if(numA<numNumerTerms)
    error('The degree of the numerator must be <= that of the denominator to use this function.')
end

%For the first pole, we have to find coefficients a of the terms
%a(1)/(x-poles(1))+a(2)/(x-poles(1))^2+...a(k(1))/(x-poles(1))^k(1)+
%The for the second pole, we find cofficients of the terms
%a(k(1)+1)/(x-poles(2))+a(k(1)+2)/(x-poles(2))^2+...a(k(1)+k(2))/(x-poles(2))^k(2)+
%and so on for all of the other poles.

%To simplify things, we precompute all polynomials (x-pole(i))^n for n=1 to
%k(i). For easy access (though not memory efficientcy, we store the
%polynomials in a 3D matrix M(:,order,i) holds the coefficients for n=order
%of pole(i). Thus, there are at most max(k)+1 coefficients.
maxK=max(k);
M=zeros(maxK+1,maxK,numPoles);
for curPole=1:numPoles
    %The highest order term is last. That is the opposite of the format in
    %Matlab's polyval function.
    polyCur=[-poles(curPole);1];
    M(1:2,1,curPole)=polyCur;
    for kCur=2:k(curPole)
        M(1:(kCur+1),kCur,curPole)=conv(polyCur,M(1:kCur,kCur-1,curPole));
    end
end 

%Having precomputed all of the powers of the individual terms, we construct
%the matrix to solve for the partial fraction coefficients. The columns
%will be multipling the a(1)...a(numA) terms. The rows will hold
%coefficients of the x terms.
C=zeros(numA,numA);
%The first element of b is the ones term.
b=[flipud(numerPoly(:));zeros(numA-numNumerTerms,1)];

%We are solving C*a=b for the a vector.
aIdx=1;
for curPole=1:numPoles
    %For all of the coefficients associated with this pole, there will be a
    %polynomial involving the maximum k for all of the other poles times
    %something involving this pole. Here, we compute the polynomial that
    %involves all of the other powers.
    polyOther=1;
    for otherPole=[1:(curPole-1),(curPole+1):numPoles]
        kOther=k(otherPole);

        polyOther=conv(polyOther,M(1:(kOther+1),kOther,otherPole));
    end
    polyOtherLen=length(polyOther);

    for kCur=1:k(curPole)
        %If the denominator in the partial fraction has order kCur, then
        %after dividing that out of the numerator, we are left with
        %k(curPole)-kCur terms.
        kLeft=k(curPole)-kCur;
        numTerms=polyOtherLen+kLeft;

        if(kLeft>0)
            C(1:numTerms,aIdx)=conv(M(1:(kLeft+1),kLeft,curPole),polyOther);
        else
            C(1:numTerms,aIdx)=polyOther;
        end
        aIdx=aIdx+1;
    end
end

coeffs=C\b;
if(nargout>1)
    polesRep=zeros(numA,1);
    idxStart=1;
    for curPole=1:numPoles
        sel=idxStart:(idxStart+k(curPole)-1);
        polesRep(sel)=poles(curPole);
        idxStart=idxStart+k(curPole);
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
