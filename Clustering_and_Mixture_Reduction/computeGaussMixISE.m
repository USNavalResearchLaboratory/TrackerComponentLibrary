function [ISEVal,PDFVals12,PDFVals22,PDFVals11,Jhr,Jrr,Jhh]=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2,outputType)
%%COMPUTEGAUSSMIXISE Given two Gaussian mixture distributions, compute the
%             integral squared error (ISE), also known as the integral
%             squared distance, (ISD), between the two distributions. For
%             distributions f_1 and f_2, the ISE is defined as
%             integral_x (f_1(x)-f_2(x))^2 dx The ISE is used as a cost
%             function for Gaussian mixture reduction in [1].
%
%INPUTS: w1 The N1X1 or 1XN1 vector of weights for the first Gaussian
%           mixture. All w1>0 and sum(w1)=1.
%       mu1 The xDimXN1 set of mean vectors for the first Gaussian mixture
%           distribution.
%        P1 The xDimXxDimXN1 set of positive definite covariance matrices
%           for the first Gaussian mixture distribution.
% w2, mu2, P2, The length N2, xDimXN2, and xDimXxDimXN2 set of weights,
%           mean vectors and positive-definite covariance matrices for the
%           second Gaussian mixture distribution.
% outputType An optional parameter indicating the type of output to provie.
%           Possible values are:
%           0 (The default if omitted or an empty matrix is passed) Return
%             the ISE in ISEVal.
%           1 Return the ISE but without the term involving only values
%             related to w1,mu1, and P1. In algorithms when optimizing over
%             w2, mu2, and P2, this term is constant.
%           2 Return the normalized ISE (NISE). The NISE is normalized such
%             that it ranges from 0 to 1, whereas the standard ISE is
%             bounded at the lower end by 0 and is unbounded at the upper
%             end.
%       
%OUTPUTS: ISEVal The scalar value of the ISE or the normalized ISE.
%      PDFVals12 A matrix such that the value in element (i,j) is
%                N(mu1(:,i);mu2(:,j),P1(:,:,i)+P2(:,:,j)), where N
%                indicates the multivariate Gaussian PDF evaluated at the
%                first argument with the second and third arguments being
%                the mean and covarince matrix. These values can be useful
%                when evaluating gradients of the ISE for optimization.
%      PDFVals22 A matrix such that the value in element (i,j) is
%                N(mu2(:,i);mu2(:,j),P2(:,:,i)+P2(:,:,j)).
%      PDFVals11 A matrix such that the value in element (i,j) is
%                N(mu1(:,i);mu1(:,j),P1(:,:,i)+P1(:,:,j)).
%  Jhr, Jrr, Jhh In Chapter 3 of [1], the ISE is expressed in terms of
%                scalar Jhr, Jrr, and Jhh terms. These are thos terms. Note
%                that if outputType=1, then Jhh is not computed and zero is
%                just returned in its place.
%
%Explicit expressions for the ISE of two Gaussian mixtures are given in
%Section 3.3.2 of [1].
%
%EXAMPLE 1:
%Here, we consider a 1D example and then the mixture missing a few
%components and also the moment-matched Gaussian approximation. In each
%instance, the ISE increases. 
% w1=[0.03,0.18,0.12,0.19,0.02,0.16,0.06,0.1,0.08,0.06];
% n1=length(w1);
% mu1=[1.45,2.20,0.67,0.48,1.49,0.91,1.01,1.42,2.77,0.89];
% P1=[0.0487,0.0305,0.1171,0.0174,0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P1=reshape(P1,[1,1,n1]);
% 
% %The second PDF is the first with the five least-weight components
% %deleted.
% w2=[0.18,0.12,0.19,0.16,0.1,0.08];
% w2=w2/sum(w2);
% n2=length(w2);
% mu2=[2.20,0.67,0.48,0.91,1.42,2.77];
% P2=[0.0305,0.1171,0.0174,0.0102,0.0380,0.0115];
% P2=reshape(P2,[1,1,n2]);
% 
% %The third PDF is the second PDF with two more components deleted.
% w3=[0.18,0.12,0.19,0.16];
% w3=w3/sum(w3);
% n3=length(w3);
% mu3=[2.20,0.67,0.48,0.917];
% P3=[0.0305,0.1171,0.0174,0.0102];
% P3=reshape(P3,[1,1,n3]);
% 
% %The fourth PDF is just the Gaussian matching the first two moments of
% %the original PDF.
% w4=1;
% [mu4,P4]=calcMixtureMoments(mu1,w1,P1);
% 
% ISEVal11=computeGaussMixISE(w1,mu1,P1,w1,mu1,P1)
% ISEVal12=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2)
% ISEVal13=computeGaussMixISE(w1,mu1,P1,w3,mu3,P3)
% ISEVal14=computeGaussMixISE(w1,mu1,P1,w4,mu4,P4)
% %With fewer elements, the ISE increases. ISEVals11 should be within
% finite precision limits of 0.
% %We can visualize the above PDFs.
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w1,mu1,P1);
% PDFVals2=GaussianMixtureD.PDF(xVals,w2,mu2,P2);
% PDFVals3=GaussianMixtureD.PDF(xVals,w3,mu3,P3);
% PDFVals4=GaussianMixtureD.PDF(xVals,w4,mu4,P4);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'-r','linewidth',2)
% plot(xVals,PDFVals3,'--g','linewidth',2)
% plot(xVals,PDFVals4,'-.b','linewidth',2)
% legend('10 Elements', '5 Elements', '3 Elements','Matched Gaussian')
%
%EXAMPLE 2:
%Here is an example with bivariare Gaussian mixtures. 
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
% ISEVal11=computeGaussMixISE(w1,mu1,P1,w1,mu1,P1)
% ISEVal12=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2)
% ISEVal13=computeGaussMixISE(w1,mu1,P1,w3,mu3,P3)
%ISEVal11 will be within finite precision limits of 0. Also, ISEVal12 will
%be larger than ISEVal13, indicating that the mixture gaussian is better
%than just dropping individual components.
%
%REFERENCES:
%[1] J. L. Williams, "Gaussian mixture reduction for tracking multiple
%    maneuvering targets in clutter," Master's thesis, Air Force Institute
%    of Technology, Mar. 2003. [Online].
%    Available: http://www.dtic.mil/srch/doc?collection=t3&id=ADA415317
%
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(outputType))
    outputType=0;
end

h=length(w1);
r=length(w2);

%Expressions for Jrr, Jhr and Jhh are given in Equation 3.25 in Section
%3.3.2 of [1].

Jhh=0;
PDFVals11=[];
if(outputType~=1)
    %If outputType==1, we don't compute Jhh and PDFVals11 is just left
    %empty.
    if(nargout>2)%If we have to return PDFVals11
        for l1=1:h
            %The following line is equivalent to 
            %w1(l1)^2*GaussianD.PDF(mu1(:,l1),mu1(:,l1),2*P1(:,:,l1));
            %because the argument of the exponent is zero.
            Jhh=Jhh+w1(l1)^2*(det(4*pi*P1(:,:,l1)))^(-1/2);
            for l2=(l1+1):h
                Jhh=Jhh+2*w1(l1)*w1(l2)*GaussianD.PDF(mu1(:,l1),mu1(:,l2),P1(:,:,l1)+P1(:,:,l2));
            end
        end
    else
        PDFVals11=zeros(h,h);
        for l1=1:h
            %The following line is equivalent to 
            %GaussianD.PDF(mu1(:,l1),mu1(:,l1),2*P1(:,:,l1));
            %because the argument of the exponent is zero.
            PDFVals11(l1,l1)=(det(4*pi*P1(:,:,l1)))^(-1/2);
            Jhh=Jhh+w1(l1)^2*PDFVals11(l1,l1);
            for l2=(l1+1):h
                PDFVals11(l1,l2)=GaussianD.PDF(mu1(:,l1),mu1(:,l2),P1(:,:,l1)+P1(:,:,l2));
                PDFVals11(l2,l1)=PDFVals11(l1,l2);
                Jhh=Jhh+2*w1(l1)*w1(l2)*PDFVals11(l1,l2);
            end
        end
    end
end

if(nargout<=1)%If PDFVals12 and PDFVals22 are not desired.
    Jhr=0;
    for l=1:h
        for n=1:r
            Jhr=Jhr+w1(l)*w2(n)*GaussianD.PDF(mu1(:,l),mu2(:,n),P1(:,:,l)+P2(:,:,n));
        end
    end
    
    Jrr=0;
    for n1=1:r
        %The following line is equivalent to 
        %Jrr=Jrr+w2(n1)^2*GaussianD.PDF(mu2(:,n1),mu2(:,n1),2*P2(:,:,n1));
        %because the argument of the exponent is zero.
        Jrr=Jrr+w2(n1)^2*(det(4*pi*P2(:,:,n1)))^(-1/2);
        for n2=(n1+1):r
            Jrr=Jrr+2*w2(n1)*w2(n2)*GaussianD.PDF(mu2(:,n1),mu2(:,n2),P2(:,:,n1)+P2(:,:,n2));
        end
    end
else%If the PDFVals12 and PDFVals22 are desired.
    PDFVals12=zeros(h,r);
    PDFVals22=zeros(r,r);

    for l=1:h
        for n=1:r
            PDFVals12(l,n)=GaussianD.PDF(mu1(:,l),mu2(:,n),P1(:,:,l)+P2(:,:,n));
        end
    end

    for n1=1:r
        %The following line is equivalent to 
        %GaussianD.PDF(mu2(:,n1),mu2(:,n1),2*P2(:,:,n1));
        %because the argument of the exponent is zero.
        PDFVals22(n1,n1)=(det(4*pi*P2(:,:,n1)))^(-1/2);
        for n2=(n1+1):r
            PDFVals22(n1,n2)=GaussianD.PDF(mu2(:,n1),mu2(:,n2),P2(:,:,n1)+P2(:,:,n2));
            PDFVals22(n2,n1)=PDFVals22(n1,n2);
        end
    end

    w1=w1(:);
    w2=w2(:);
    
    Jhr=sum(vec((w1*w2').*PDFVals12));
    Jrr=sum(vec((w2*w2').*PDFVals22));
end

%Equation 3.20 in [1]. The ISE
ISEVal=Jhh+Jrr-2*Jhr;
switch(outputType)
    case 0
    case 1
    case 2%Normalize to get the NISE.
        %The NISE comes from Equation 3.26 in [1].
        ISEVal=ISEVal/(Jrr+Jhh);
    otherwise
        error('Unknown output type specified.')
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
