function [x,P,gammaVec,Gamma]=ellipsoidIntersect(xi,Pi,xj,Pj)
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
%        Gamma The mutual covariance matrix used in the fusion.
%
%The algorithm in [1] is used, (which is an expanded version of [2]), where
%the mutual mean and covariance matrix are defined.
%
%EXAMPLE:
%Here, we recreate the example in Fig. 5 of [2].
% xi=[1;-2];
% Pi=[3,0;
%     0,0.4];
% xj=[-2;-1];
% Pj=[2,-0.8;
%     -0.8,1];
% [x,P,gammaVec,Gamma]=ellipsoidIntersect(xi,Pi,xj,Pj);
% figure(1)
% clf
% hold on
% threshold=1;
% drawEllipse(xi,inv(Pi),threshold,'--r')
% drawEllipse(xj,inv(Pj),threshold,'--g')
% drawEllipse(x,inv(P),threshold,'-b')%The fused estimate.
% drawEllipse(gammaVec,inv(Gamma),threshold,'-.k')%The mutual estimate.
%
%REFERENCES:
%[1] J. Sijs and M. Lazar, "State fusion with unknown correlation:
%    Ellipsoidal intersection," Automatica, vol. 48, no. 8, pp. 1874-1878,
%    Aug. 2012.
%[2] J. Sijs, M. Lazar, and P. Bosch, "State fusion with unknown
%    correlation: Ellipsoidal intersection," in Proceedings of the 2010
%    American Control Conference, Baltimore, MD, 30 Jun. - 2 Jul. 2010.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation numbers are taken from [1].
xDim=size(Pi,1);

%Equation 12a
[Si,Di]=eig(Pi);
di=diag(Di);
diRoot=sqrt(di);
DiRoot=diag(diRoot);
DiRootInv=diag(1./diRoot);
[Sj,Dj]=eig(DiRootInv*Si\Pj*Si*DiRootInv);
dj=diag(Dj);

%Equation 12b
T=Si*DiRoot*Sj;

%Equation 10 defines \hat{P}_i as inv(T)*Pi*inv(T)' and \hat{P}_j the same
%way substituting Pj. The result is that \hat{P}_i is the identity matrix
%and \hat{P}_j is the diagonal matrix whose main diagonal is dj.

mui=T\xi;
muj=T\xj;

%DGamma in Equation 11 and for finding muGamma in Equation 14a.
dGamma=zeros(xDim,1);
for curDim=1:xDim
    if(dj(curDim)==1)
        dGamma(curDim)=1;
    elseif(dj(curDim)<1)
        dGamma(curDim)=1;
    elseif(dj(curDim)>1)
        dGamma(curDim)=dj(curDim);
    end
end
DGamma=diag(dGamma);

%Theorem 6
Gamma=T*DGamma*T';

%Equation 16b, without the eta term.
Wi=diag(1./dj-1./dGamma);
Wj=diag(1-1./dGamma);

%Make sure that zeta is large enough to make a difference.
zeta=max(max(eps(Wi(:))),max(eps(Wj(:))));
%Equation 16c.
if(abs(dj-1)>=10*zeta)
    eta=0; 
else
    eta=10*zeta;
end

%Add the eta term to Equation 16b.
Wi=Wi+eta*eye(xDim,xDim);
Wj=Wj+eta*eye(xDim,xDim);

%Equation 16a
muGamma=(Wi+Wj)\(Wi*mui+Wj*muj);
%Described after Equation 16c.
gammaVec=T*muGamma;

%Equation 8. We using pinv for better numerical stability.
P=pinv(inv(Pi)+inv(Pj)-inv(Gamma));
x=P*(Pi\xi+Pj\xj-Gamma\gammaVec);

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
