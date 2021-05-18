function distVal=GaussianSimilarity(mu1,Sigma1,mu2,Sigma2)
%%GAUSSIANSIMILARITY This is the normalized similarity measure used between
%               two Gaussian distributions in the Appendix of [1]. It
%               ranges from 0 (completely different) to 1 (same
%               distribution).
%
%INPUTS: mu1 The xDimX1 mean of the first Gaussian distribution.
%     Sigma1 The xDimXxDim covariance matrix of the first Gaussian
%            distribution.
%        mu2 The xDimX1 mean of the second Gaussian distribution.
%     Sigma2 The xDimXxDim covariance matrix of the second Gaussian
%            distribution.
%
%OUTPUTS: distVal The scalar distance between the two distributions.
%
%EXAMPLE:
%We show the similairty of identical distirbutions, extremely different
%distributions, and distributions that just differ a little bit.
% mu1=[12;3];
% mu2=[2000;8000];
% Sigma1=[10, 6;
%          6, 9];
% Sigma2=[8, -3;
%        -3, 14];
% %Identical Distributions (Similarity is 1).
% GaussianSimilarity(mu1,Sigma1,mu1,Sigma1)
% %Distributions that are very different (Similarity is 0).
% GaussianSimilarity(mu1,Sigma1,mu2,Sigma2)
% %Distributions that only differ a little bit (Similarity is about
% %0.2163).
% GaussianSimilarity(mu1,Sigma1,mu1*0.5,Sigma1*0.5)
%
%REFERENCES:
%[1] F. Faubel, J. McDonough, and D. Klakow, "The split and merge unscented
%    Gaussian mixture filter," IEEE Signal Processing Letters, vol. 16, no.
%    9, pp. 786-789, Sep. 2009.
%
%January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 16.
Sigma12Inv=inv(Sigma1)+inv(Sigma2);
%Equation 15.
mu12=Sigma12Inv\(Sigma1\mu1+Sigma2\mu2);

numDim=size(mu1,1);
z=zeros(numDim,1);

p1=GaussianD.PDF(z,mu1,Sigma1/2);
p2=GaussianD.PDF(z,mu2,Sigma2/2);
p12=GaussianD.PDFI(z,mu12,Sigma12Inv);

if(p1==0||p2==0||p12==0||p1*p2==0)
    distVal=0;
else
    distVal=sqrt(p1*p2)/p12;
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
