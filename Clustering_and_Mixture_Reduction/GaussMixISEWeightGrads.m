function [gradISE,gradJhr,gradJrr]=GaussMixISEWeightGrads(w1,w2,sqrtw2,PDFVals12,PDFVals22)
%%GAUSSMIXISEWEIGHTGRADS Compute the gradient of the non-normalized
%               integrated squared error (ISE) between two Gaussian mixture
%               PDFs with respect to the square root of the weight terms of
%               the second Gaussian mixture. This gradient can arise when
%               optimizing over the ISE with respect to the square root of
%               the weights. The square root of the weights is typically
%               used so as to ensure that no weight ever becomes negative.
%
%INPUTS: w1 The 1Xn1 or n1X1 set of weights of the first Gaussian mixture
%           distribution. It is required that all w1>=0 and sum(w1)=1.
%        w2 The 1Xn2 or n2X1 set of weights of the second Gaussian mixture
%           distribution. It is required that all w2>=0 and sum(w2)=1.
%    sqrtw2 This is just sqrt(w2) and is the term with which the gradient
%           is being taken.
% PDFVals12 A matrix such that the value in element (i,j) is
%           N(mu1(:,i);mu2(:,j),P1(:,:,i)+P2(:,:,j)), where N indicates the
%           multivariate Gaussian PDF evaluated at the first argument with
%           the second and third arguments being the mean and covarince
%           matrix. mu1 and P1 are the means and covariance matrices of the
%           first Gaussian mixture distribution and mu2 and P2 are the same
%           for the second Gaussian mixture distribution. This parameter is
%           returned by computeGaussMixISE.
% PDFVals22 A matrix such that the value in element (i,j) is
%           N(mu2(:,i);mu2(:,j),P2(:,:,i)+P2(:,:,j)). This parameter is
%           returned by computeGaussMixISE.
%
%OUTPUTS: gradISE The 1Xn2 set of derivatives of the ISE wtih respect to
%                 the elements of q (the square roots of w2).
% gradJhr, gradJrr In Chapter 3 of [1], the ISE is expressed in terms of
%                 Jhr and Jrr terms. These are the 1Xn2 gradients of those
%                 terms.
%
%Formule for gradJhr and gradJrr are Equation 3.31 in Section 3.3.3.1 of
%[1]. They relate to the ISE via Equation 3.20. See the function
%computeGaussMixISE to compute the ISE.
%
%EXAMPLE:
%In this example with a scalar PDF, we verify that the gradient obtained
%from this function is consistent with numerical differentiation.
% w1=[0.03,0.18,0.12,0.19,0.02,0.16,0.06,0.1,0.08,0.06];
% n1=length(w1);
% mu1=[1.45,2.20,0.67,0.48,1.49,0.91,1.01,1.42,2.77,0.89];
% P1=[0.0487,0.0305,0.1171,0.0174,0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P1=reshape(P1,[1,1,n1]);
% 
% %The second PDF is the first with the five least-weight components deleted.
% w2=[0.18,0.12,0.19,0.16,0.1,0.08];
% w2=w2/sum(w2);
% n2=length(w2);
% mu2=[2.20,0.67,0.48,0.91,1.42,2.77];
% P2=[0.0305,0.1171,0.0174,0.0102,0.0380,0.0115];
% P2=reshape(P2,[1,1,n2]);
% 
% [ISEVal,PDFVals12,PDFVals22]=computeGaussMixISE(w1,mu1,P1,w2,mu2,P2);
% q2=sqrt(w2);
% epsVal=1e-9;
% gradISENum=zeros(1,n2);
% for k=1:n2
%     q2Cur=q2;
%     q2Cur(k)=q2Cur(k)+epsVal;
%     ISEValCur=computeGaussMixISE(w1,mu1,P1,q2Cur.^2,mu2,P2);
%     gradISENum(k)=(ISEValCur-ISEVal)/epsVal;
% end
% gradISE=GaussMixISEWeightGrads(w1,w2,PDFVals12,PDFVals22);
% RelErr=max(abs((gradISENum-gradISE)./gradISENum))
%The relative error will be about 2.1628e-6, which indicates good numeric
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

q=sqrtw2;

%Both formulae are Equation 3.31 in Section 3.3.3.1 of [1].
gradJhr=2*q(:).'.*sum(bsxfun(@times,w1(:),PDFVals12),1);
gradJrr=4*q(:).'.*sum(bsxfun(@times,w2(:),PDFVals22),1);

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
