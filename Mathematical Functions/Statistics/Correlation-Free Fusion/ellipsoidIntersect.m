function [x,P,gammaVec,GammaMat]=ellipsoidIntersect(xi,Pi,xj,Pj)
%%ELLIPSOIDINTERSECT Perform ellipsoidal intersection to merge two Gaussian
%                    estimates having an unknown correction. Compare this
%                    function to covarianceIntersect and
%                    covarianceIntersectSCI.
%
%INPUTS: xi,Pi The nX1 mean vector and nXn positive definite covariance
%              matrix of the first estimate.
%        xj,Pj The nX1 mean vector and nXn positive definite covariance
%              matrix of the second estimate.
%
%OUTPUTS: x, P The nX1 mean and nXn covariance matrix of the estimate fused
%              using ellipsoidal intersection.
%     gammaVec The mutual mean vector used in the fusion.
%     GammaMat The mutual covariance matrix used in the fusion.
%
%The algorithm in [1] is used, where the mutual mean and covariance matrix
%are defined.. However, rather than using diagonal loading to deal with the
%inverse of a possibly singular matrix in Equation 13b, we use a
%pseudoinverse.
%
%EXAMPLE:
%Here, we use the numbers from the example in [1] and we plot ellipsoids.
%Fig. 5 in the paper appears to be in error.
%is in error.
% xi=[1;-2];
% Pi=[3,0;0,0.4];
% xj=[-2;-1];
% Pj=[2,-0.8;-0.8,1];
% 
% [x,P]=ellipsoidIntersect(xi,Pi,xj,Pj);
% figure()
% clf
% hold on
% drawEllipse(xi,inv(Pi),[],'--r')
% drawEllipse(xj,inv(Pj),[],'--g')
% drawEllipse(x,inv(P),[],'-b')
%
%REFERENCES:
% [1] J. Sijs, M. Lazar, and P. Bosch, "State fusion with unknown
%     correlation: Ellipsoidal intersection," in Proceedings of the 2010
%     American Control Conference, Baltimore, MD, 30 Jun. - 2 Jul. 2010.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Si,Di]=eig(Pi);
[Sj,Dj]=eig(Pj);

PiInv=inv(Pi);
PjInv=inv(Pj);

%Equation 11b
DGamma=diag(max(diag(Dj),1));
%Equation 11a, mutual covariance
GammaMat=Si*sqrt(Di)*(Sj*DGamma/Sj)*sqrt(Di)/Si;
GammaMatInv=inv(GammaMat);

%Before Equation 14. 
Wi=PiInv-GammaMatInv;
Wj=PjInv-GammaMatInv;

%We choose to use a pseudoinverse rather than diagonal loading to implement
%Equation 13b.
gammaVec=inv(Wi+Wj)*(Wi*xi+Wj*xj);

%Equation 8. We using pinv for better numerical stability.
P=pinv(PiInv+PjInv-GammaMatInv);
x=P*(PiInv*xi+PjInv*xj-GammaMatInv*gammaVec);

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
