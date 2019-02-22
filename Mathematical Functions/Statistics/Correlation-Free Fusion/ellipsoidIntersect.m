function [P,Gamma,x,gammaVec]=ellipsoidIntersect(Pi,Pj,xi,xj,zeta)
%%ELLIPSOIDINTERSECT Perform ellipsoidal intersection to merge two Gaussian
%                    estimates having an unknown correlation. Compare this
%                    function to covarianceIntersect and
%                    covarianceIntersectSCI.
%
%INPUTS: Pi, Pj The nXn positive definite covariance matrices of the first
%               and second estimates.
%        xi, xj The nX1 mean vectors of the first and second estimates. If
%               one does not want x or gammaVec on the output, then these
%               can be omitted.
%          zeta An optional parameter that is a small value that is used to
%               determine approximate equality of some eigenvalues. If
%               omitted or an empty matrix is passed, the smallest value
%               possible that can make a difference is chosen.
%
%OUTPUTS: P The nXn covariance matrix of the fused estimate.
%     Gamma The mutual covariance matrix used in the fusion.
%         x The nX1 mean value of the fused estimate.
%  gammaVec The mutual mean vector used in the fusion.
%
%The algorithm in [1] is used, (which is an expanded discussion of what is
%in [2]), where the mutual mean and covariance matrix are defined.
%
%EXAMPLE 1:
%Here, we recreate the example in Fig. 5 of [2].
% xi=[1;-2];
% Pi=[3,0;
%     0,0.4];
% xj=[-2;-1];
% Pj=[2,-0.8;
%     -0.8,1];
% [P,Gamma,x,gammaVec]=ellipsoidIntersect(Pi,Pj,xi,xj);
% figure(1)
% clf
% hold on
% threshold=1;
% drawEllipse(xi,inv(Pi),threshold,'--r','linewidth',2)
% drawEllipse(xj,inv(Pj),threshold,'--g','linewidth',2)
% drawEllipse(x,inv(P),threshold,'-b','linewidth',2)%The fused estimate.
% drawEllipse(gammaVec,inv(Gamma),threshold,'-.k','linewidth',2)%The mutual estimate
% axis([-3.5 4.5 -3 1])
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%legend({'\epsilon_{x_i,P_i}(x)','\epsilon_{x_j,P_j}(x)','\epsilon_{x_\gamma,\Gamma}(x)','\epsilon_{x_{if},P_{if}}(x)'},'FontSize',18)
%
%EXAMPLE 2:
%Here, we recreate Fig. 2 of [1] without the covariance intersection
%ellipse.
% xi=[0.5;1];
% Pi=[2.5,-1;
%     -1,1.2];
% xj=[2;1];
% Pj=[0.8,-0.5;
%     -0.5,4];
% zeta=0.1;
% [P1,~,x1]=ellipsoidIntersect(Pi,Pj,xi,xj,zeta);
% zeta=1e-6;
% [P2,~,x2]=ellipsoidIntersect(Pi,Pj,xi,xj,zeta);
% figure(1)
% clf
% hold on
% threshold=0.9;
% drawEllipse(xi,inv(Pi),threshold,'--r','linewidth',2)
% drawEllipse(xj,inv(Pj),threshold,'--g','linewidth',2)
% %The fused estimates.
% drawEllipse(x1,inv(P1),threshold,'-b','linewidth',2)
% drawEllipse(x2,inv(P2),threshold,'-b','linewidth',2)
% axis([-1.5 3 -1 3.1])
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend({'\epsilon_{x_i,P_i}(x)','\epsilon_{x_j,P_j}(x)','\epsilon_{x_f,P_{f}}(x)'},'FontSize',18,'Location','SouthWest')
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

if(nargin<5||isempty(zeta))
    %Choose the smallest zeta that is large enough to make a difference.
    zeta=max(max(eps(Wi(:))),max(eps(Wj(:))));
end
%Equation 16c.
if(abs(dj-1)>=10*zeta)
    eta=0; 
else
    eta=zeta;
end

%Add the eta term to Equation 16b.
Wi=Wi+eta*eye(xDim,xDim);
Wj=Wj+eta*eye(xDim,xDim);

%Equation 8. We using pinv for better numerical stability.
P=pinv(inv(Pi)+inv(Pj)-inv(Gamma));

if(nargout>2)
    %Equation 13
    mui=T\xi;
    muj=T\xj;
    %Equation 16a
    muGamma=(Wi+Wj)\(Wi*mui+Wj*muj);
    %Described after Equation 16c.
    gammaVec=T*muGamma;

    x=P*(Pi\xi+Pj\xj-Gamma\gammaVec);
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
