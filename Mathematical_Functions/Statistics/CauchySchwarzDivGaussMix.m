function DCS=CauchySchwarzDivGaussMix(w1,mu1,P1,w2,mu2,P2)
%%CAUCHYSCHWARZDIVGAUSSMIX Compute the Cauchy-Scharz divergence between two
%               Gaussian mixture PDFs. Unlike the Kullback-Leibler (KL)
%               divergence, there is a closed-form expression for the
%               Cauchy-Schwarz divergence between Gaussian mixtures. Unlike
%               the KL divergence, the Cauchy-Schwarz divergence is
%               commutative: reversing the order of the mixtures as inputs
%               does not change the result.
%
%INPUTS: w1 The 1Xn1 or n1X1 set of weights of the first Gaussian mixture's
%           components such that all w>=0 and sum(w)=1.
%       mu1 The dXn1 set of mean values of the first Gaussian mixture's
%           components.
%        P1 The dXxdXn1 set of covariance matrices of the first Gaussian
%           mixture's components.
% w2,mu2,P2 The equivalent of w1, mu1 and P1 for the second mixtures
%           components. There length-n2, dXn2, and dXdXn2 in size.
%
%OUTPUTS: DCS The Cauchy-Scharz divergence between the two gaussian
%             mixtures.
%
%This function implements Equation 3 of [1].
%
%EXAMPLE:
%IN this example, we consider compare a 3-component Gaussian mixture
%against itself (which gives a Cauchy-Schwarz divergence of zero), as well
%as against a mixture missing one component and against a moment-matched
%single Gaussian PDF. with fewer Gaussians, the value of the Cauchy-Schwarz
%divergence increases.
% w1=[0.25;0.5;0.25];
% mu1=zeros(2,2);
% mu1(:,1)=[1;-1];
% mu1(:,2)=[-1;1];
% mu1(:,3)=[0;0];
% P1=zeros(2,2,2);
% P1(:,:,1)=[4/9,  14/45;
%            14/45,4/9];
% P1(:,:,2)=[4/9, 0;
%            0, 4/9];
% P1(:,:,3)=[2/9, -1/9;
%           -1/9, 3/9];
% 
% %The second distribution just throws out the first component.
% w2=w1(2:3);
% w2=w2/sum(w2);
% mu2=mu1(:,2:3);
% P2=P1(:,:,2:3);
% 
% %The third distribution replaces the original distribution with a single
% %moment-matched Gaussian.
% w3=1;
% [mu3,P3]=calcMixtureMoments(mu1,w1,P1);
% div11=CauchySchwarzDivGaussMix(w1,mu1,P1,w1,mu1,P1)
% div12=CauchySchwarzDivGaussMix(w1,mu1,P1,w2,mu2,P2)
% div13=CauchySchwarzDivGaussMix(w1,mu1,P1,w3,mu3,P3)
%
%REFERENCES:
%[1] K. Kampa, E. Hasanbelliu, and J. C. Principe, "Closed-form Cauchy-
%    Schwarz PDF divergence for mixture of Gaussians," in Proceedings of
%    the International Joint Conference on Neural Networks, San Jose, CA,
%    31 Jul.-5 Aug. 2011, pp. 2578-2585.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n1=length(w1);
n2=length(w2);

%zmk, zmm and zkk are given after Equation 3.
zmk=zeros(n1,n2);
zmm=zeros(n1,n1);
zkk=zeros(n2,n2);

for m=1:n1
    for k=1:n2
        zmk(m,k)=GaussianD.PDF(mu1(:,m),mu2(:,k),P1(:,:,m)+P2(:,:,k));
    end
end

%We only need to fill in the lower-triangular parts of zmm and zkk
for m1=1:n1
    for m2=1:(m1-1)
        zmm(m1,m2)=GaussianD.PDF(mu1(:,m1),mu1(:,m2),P1(:,:,m1)+P1(:,:,m2));
    end
end

for k1=1:n2
    for k2=1:(k1-1)
        zkk(k1,k2)=GaussianD.PDF(mu2(:,k1),mu2(:,k2),P2(:,:,k1)+P2(:,:,k2));
    end
end

term1=0;
for m=1:n1
    for k=1:n2
        term1=term1+w1(m)*w2(k)*zmk(m,k);
    end
end

term2=0;
for m=1:n1
    term2=term2+w1(m)^2/sqrt(det(2*pi*P1(:,:,m)));
end

term3=0;
for m1=1:n1
    for m2=1:(m1-1)
        term3=term3+w1(m1)*w1(m2)*zmm(m1,m2);
    end
end

term4=0;
for k=1:n2
    term4=term4+w2(k)^2/sqrt(det(2*pi*P2(:,:,k)));
end

term5=0;
for k1=1:n2
    for k2=1:(k1-1)
        term5=term5+w2(k1)*w2(k2)*zkk(k1,k2);
    end
end

%Equation 3 in [1].
DCS=-log(term1)+(1/2)*log(term2+2*term3)+(1/2)*log(term4+2*term5);

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
