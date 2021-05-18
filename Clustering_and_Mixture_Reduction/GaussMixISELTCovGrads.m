function [gradISE,gradJhr,gradJrr]=GaussMixISELTCovGrads(w1,mu1,P1,w2,mu2,P2,L2,PDFVals12,PDFVals22)
%%GAUSSMIXISELTCOVGRADS Compute the gradient of the non-normalized
%               integrated squared error (ISE) between two Gaussian mixture
%               PDFs with respect to the elements of square-roots of the
%               covariance matrices of the individual components of the
%               second Gaussian mixture. For a covariance matrix P, the
%               square root L is such that P=L*L'. However, L does not need
%               to be lower-triangular. This gradient can arise when
%               optimizing over the ISE with respect to the covariance
%               matrices. Optimizing over the square-root matrices rather
%               than the matrices themselves ensures that the covariance
%               matrices are always positive (semi)definite.
%
%INPUTS: w1 The N1X1 or 1XN1 vector of weights for the first Gaussian
%           mixture. All w1>0 and sum(w1)=1.
%       mu1 The xDimXN1 set of mean vectors for the first Gaussian mixture
%           distribution.
%        P1 The xDimXxDimXN1 set of positive definite covariance matrices
%           for the first Gaussian mixture distirbution.
% w2, mu2, P2, The length N2, xDimXN2, and xDimXxDimXN2 set of weights,
%           mean vectors and positive-definite covariance matrices for the
%           second Gaussian mixture distribution.
%        L2 The xDimXxDimXN2 set of square roots of the covariance matrices
%           in P2 such that P2(:,:,i)=L2(:,:,i)*L2(:,:,i)'.
% PDFVals12 A matrix such that the value in element (i,j) is
%           N(mu1(:,i);mu2(:,j),P1(:,:,i)+P2(:,:,j)), where N indicates the
%           multivariate Gaussian PDF evaluated at the first argument with
%           the second and third arguments being the mean and covarince
%           matrix. This parameter is returned by computeGaussMixISE.
% PDFVals22 A matrix such that the value in element (i,j) is
%           N(mu2(:,i);mu2(:,j),P2(:,:,i)+P2(:,:,j)). This parameter is
%           returned by computeGaussMixISE.
%
%OUTPUTS: gradISE The xDimXxDimXN2 set of derivatives of the ISE with
%                 respect to the elements of L2 (the square roots of P2).
% gradJhr, gradJrr In Chapter 3 of [1], the ISE is expressed in terms of
%                 Jhr and Jrr terms. These are the xDimXxDimXN2 gradients
%                 of those terms.
%
%Formule for gradJhr and gradJrr are Equation 3.45 in Section 3.3.3.3 of
%[1]. They relate to the ISE via Equation 3.20. See the function
%computeGaussMixISE to compute the ISE.
%
%EXAMPLE 1:
%In this example with a scalar PDF, we verify that the gradient obtained
%from this function is consistent with numerical differentiation.
% w1=[0.03,0.18,0.12,0.19,0.02,0.16,0.06,0.1,0.08,0.06];
% n1=length(w1);
% mu1=[1.45,2.20,0.67,0.48,1.49,0.91,1.01,1.42,2.77,0.89];
% P1=[0.0487,0.0305,0.1171,0.0174,0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P1=reshape(P1,[1,1,n1]);

% %The second PDF is the first with the five least-weight components deleted.
% w2=[0.18,0.12,0.19,0.16,0.1,0.08];
% w2=w2/sum(w2);
% n2=length(w2);
% mu2=[2.20,0.67,0.48,0.91,1.42,2.77];
% P2=[0.0305,0.1171,0.0174,0.0102,0.0380,0.0115];
% P2=reshape(P2,[1,1,n2]);
% 
% [ISEVal,PDFVals12,PDFVals22]=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2);
% L2=sqrt(P2);
% epsVal=1e-8;
% gradISENum=zeros(1,1,n2);
% for k=1:n2
%     L2Cur=L2;
%     L2Cur(1,1,k)=L2Cur(1,1,k)+epsVal;
%     ISEValCur=computeGaussMixISE(w1,mu1,P1,w2,mu2,L2Cur.^2);
%     gradISENum(k)=(ISEValCur-ISEVal)/epsVal;
% end
% gradISE=GaussMixISELTCovGrads(w1,mu1,P1,w2,mu2,P2,L2,PDFVals12,PDFVals22);
% RelErr=max(abs((gradISENum(:)-gradISE(:))./gradISENum(:)))
%The relative error will be about 7.162e-6, which indicates good numeric
%agreement.
%
%EXAMPLE 2:
%In this example with a bivariate PDF, we verify that the gradient obtained
%from this function is consistent with numerical differentiation.
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
% n2=length(w2);
% [ISEVal,PDFVals12,PDFVals22]=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2);
% 
% numDim=size(mu2,1);
% L2=zeros(numDim,numDim,n2);
% for j=1:n2
%     L2(:,:,j)=chol(P2(:,:,j),'lower'); 
% end
% 
% epsVal=1e-8;
% gradISENumDiff=zeros(numDim,numDim,n2);
% for j=1:n2
%     for curEl=1:(numDim^2)
%         [i1,i2]=ind2sub([numDim,numDim],curEl);
% 
%         LCur=L2(:,:,j);
%         LCur(i1,i2)=LCur(i1,i2)+epsVal;
%         
%         PCur=P2;
%         PCur(:,:,j)=LCur*LCur';
%         ISEValCur=computeGaussMixISE(w1,mu1,P1,w2,mu2,PCur);
%         gradISENumDiff(i1,i2,j)=(ISEValCur-ISEVal)/epsVal;
%     end
% end
% gradISE=GaussMixISELTCovGrads(w1,mu1,P1,w2,mu2,P2,L2,PDFVals12,PDFVals22);
% RelErr=max(max(abs((gradISENumDiff-gradISE)./gradISENumDiff)))
%The relative error will be about 6.99e-7, which indicates good numeric
%agreement.
%
%REFERENCES:
%[1] J. L. Williams, "Gaussian mixture reduction for tracking multiple
%    maneuvering targets in clutter," Master's thesis, Air Force Institute
%    of Technology, Mar. 2003. [Online].
%    Available: http://www.dtic.mil/srch/doc?collection=t3&id=ADA415317
%
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

Nh=length(w1);
Nr=length(w2);
xDim=size(mu1,1);

%If only the covariance matrices are given.
if(nargin<7||isempty(L2))
    L2=zeros(xDim,xDim,Nr);
    for i=1:Nr
        L2(:,:,i)=chol(P2(:,:,i),'lower');
    end
end

%If only the lower-triangular Cholesky decompositions of the covariance
%matrices are given.
if(isempty(P2))
    P2=zeros(xDim,xDim,Nr);
    for i=1:Nr
        P2(:,:,i)=L2(:,:,i)*L2(:,:,i)';
    end
end

%The formulae for the gradient of Jhr and Jrr are given in Equation 3.45 of
%Section 3.3.3.3 of [1].
gradJhr=zeros(xDim,xDim,Nr);
gradJrr=zeros(xDim,xDim,Nr);

for j=1:Nr
    for i=1:Nh
        diff=mu1(:,i)-mu2(:,j);
        PSum=P1(:,:,i)+P2(:,:,j);
        PSumInv=inv(PSum);
        
        gradJhr(:,:,j)=gradJhr(:,:,j)+w1(i)*PDFVals12(i,j)*PSumInv*(diff*diff'-PSum)*PSumInv*L2(:,:,j);
    end
    gradJhr(:,:,j)=w2(j)*gradJhr(:,:,j);
end

for j=1:Nr
    for i=1:Nr
        diff=mu2(:,i)-mu2(:,j);
        PSum=P2(:,:,i)+P2(:,:,j);
        PSumInv=inv(PSum);
        
        gradJrr(:,:,j)=gradJrr(:,:,j)+w2(i)*PDFVals22(i,j)*PSumInv*(diff*diff'-PSum)*PSumInv*L2(:,:,j);
    end
    gradJrr(:,:,j)=2*w2(j)*gradJrr(:,:,j);
end

gradISE=gradJrr-2*gradJhr;
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
