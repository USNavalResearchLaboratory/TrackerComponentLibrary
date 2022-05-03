function FIM=observedFisherInfo(z,RInv,h,JacobMat,HessMat)
%%OBSERVEDFISHERINFO Assuming that a linear or nonlinear measurement is
%           corrupted with zero-mean Gaussian noise, the observed Fisher
%           information matrix (FIM) has a standard form in terms of the
%           values of the measurement function and its first and second
%           derivatives. This function takes those values and returns the
%           observed FIM. Summing the FIMs from multiple simultaneous
%           independent measurements or measurement components returns the
%           observed FIM for the fused measurement. The inverse of the FIM
%           is the Cramér-Rao lower bound (CRLB). If only a single
%           measurement is considered, and h=z, then h, z, and HessMat can
%           all be omitted. Usualy, there is no benefit to including the
%           terms.
%
%INPUTS: z The zDimX1 measurement, or if multiple measurements of the same
%          time are to the zDimXnumMeas matrix of those measurements. This
%          can be omitted if one just wants the FIM without this.
%     RInv The zDimXzDim inverse of the covariance matrix associated with
%          the multivariate Gaussian noise corrupting z, or if multiple
%          measurements are to be fused AND RInv differs among them, then a
%          zDimXzDimXnumMeas collection of all of the inverse matrices. If
%          z is omitted and multiple measurements are fused, then RInv MUST
%          be specified as a zDimXzDimXnumMeas matrix, not as a single
%          zDimXzDim matrix.
%        h The zDimX1 value of the xDimX1 state converted into the
%          measurement domain. If z is omitted, then this is not needed.
% JacobMat The zDimXxDim Jacobian matrix of derivatives of the measurement
%          function h taken with respect to the elements of the target
%          state. This is assumed the same for all measurements fused by
%          this function.
%  HessMat The xDimXxDimXzDim matrix of second derivatives of the
%          measurement function h with respect to the elements of the state
%          x. HessMat(i,j,k) is the Hessian for the kth measurement
%          component with derivatives taken with respect to elements i and
%          j of the x vector. i and j can be equal. Note that all
%          HessMat(:,:,k) are symmetric matrices. In 3D, the order of the
%          second derivatives in each submatrix is of the form:
%                  [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
%                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
%                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
%
%OUTPUTS: FIM The xDimXxDim observed Fisher information matrix.
%
%The FIM and CRLB and in many statistics texts. When considering target
%tracking, one can look at Chapter 2.7.2 of [1]. Since no expectation is
%taken for the observed Fisher information matrix, only the form in terms
%of second derivatives in [1] can be used, not the form in terms of an
%outer product of first derivatives. The use of the inverse of the observed
%FIM in characterizing the accuracy of ML estimates is discussed in [2].
%
%The observed FIM is simply the negative of the matrix of second
%derivatives of the logarithm of the likelihood function. In the problem at
%hand:
%nabla_{x}(\nabla_x)'log(p(z|x))
%where here p(z|x)=1/sqrt(det(2*pi*R))*exp(-(1/2)*(z-h(x))'*inv(R)*(z-h(x))
%The gradient of the logarithm of the likelihood function is
%nabla_{x}log(p(z|x))=-H'*inv(R)(h(x)-z)
%where H=nabla_{x} h(x)'
%The matrix of second derivatives is thus
%nabla_{x}(\nabla_x)'log(p(z|x))=-H'*inv(R)*H-C
%where the jth column of C is given by
%C(:,j)=((\partial / \partial x_j)H')*inv(R)*(h(x)-z)
%
%EXAMPLE 1:
%In this example, we consider how well the observed FIM can be used as the
%covariance matrix of a fused measurement. In this instance, with all of
%the measurement being the same accuracy, we just average the measurements
%to get a fused measurement. We then find the NEES both with and without
%using the Hessian term. It is seen that when using the Hessian term, the
%NEES is closet to 1 than when not using the Hessian term.
% numMCRuns=10000;
% numMeas=3;
% sigmaR=100;
% sigmaAz=3*(pi/180);
% sigmaEl=3*(pi/180);
% SR=diag([sigmaR;sigmaAz;sigmaEl]);
% R=SR*SR';
% RInv=inv(R);
% 
% zTrue=[30e3;60*(pi/180);3*(pi/180)];
% systemType=0;
% useHalfRange=true;
% xTrue=spher2Cart(zTrue,systemType,useHalfRange);
% 
% NEESWithHess=0;
% NEESWithoutHess=0;
% for k=1:numMCRuns
%     zMeas=zeros(3,numMeas);
%     for curMeas=1:numMeas
%         zMeas(:,curMeas)=zTrue+SR*randn(3,1);
%     end
%      xAvg=mean(spher2Cart(zMeas,systemType,useHalfRange),2);
%      
%      h=Cart2Sphere(xAvg,systemType,useHalfRange);
%      JacobMat=calcSpherJacob(xAvg,systemType,useHalfRange);
%      HessMat=calcSpherHessian(xAvg,systemType,useHalfRange);
% 
%      FIMWithHess=observedFisherInfo(zMeas,RInv,h,JacobMat,HessMat);
%      FIMWithoutHess=numMeas*observedFisherInfo([],RInv,[],JacobMat,HessMat);
%      diff=xAvg-xTrue;
%      NEESWithHess=NEESWithHess+diff'*FIMWithHess*diff;
%      NEESWithoutHess=NEESWithoutHess+diff'*FIMWithoutHess*diff;
% end
% NEESWithHess=NEESWithHess/(3*numMCRuns)
% NEESWithoutHess=NEESWithoutHess/(3*numMCRuns)
%
%EXAMPLE 2:
%In this instance, we used the observed Fisher information as a covariance
%matrix of a single spherical measurement in the absence of knowing the
%truth. Thus, we just take h(z)=z and can omit the matrix of second
%derivatives. Evaluating the NEES, one sees it is consistent (near 1, maybe
%like 0.99 or 1.001). Of course, at higher angular noise levels, a debiased
%function like spher2CartTaylor can perform better.
% numMCRuns=10000;
% sigmaR=10;
% sigmaAz=0.1*(pi/180);
% sigmaEl=0.1*(pi/180);
% SR=diag([sigmaR;sigmaAz;sigmaEl]);
% R=SR*SR';
% RInv=inv(R);
% 
% zTrue=[1e3;60*(pi/180);3*(pi/180)];
% systemType=0;
% useHalfRange=true;
% xTrue=spher2Cart(zTrue,systemType,useHalfRange);
% 
% NEES=0;
% for k=1:numMCRuns
%     zMeas=zTrue+SR*randn(3,1);
%     xConv=spher2Cart(zMeas,systemType,useHalfRange);
%     JacobMat=calcSpherJacob(xConv,systemType,useHalfRange);
%     
%     invCRLB=observedFisherInfo([],RInv,[],JacobMat);
%     
%     diff=xConv-xTrue;
%     NEES=NEES+diff'*invCRLB*diff;
% end
% NEES=NEES/(3*numMCRuns)
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation: Theory, Algorithms and
%    Software. New York: John Wiley and Sons, 2001.
%[2] B. Efron and D. Hinkley, "Assessing the accuracy of the maximum
%    likelihood estimator: Observed versus expected Fisher information,"
%    Department of Statistics, Stanford, University, Tech. Rep. 108, 8 Mar.
%    1978.
%
%December 2020 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isempty(z))
    numZ=size(z,2);
    
    xDim=size(HessMat,1);
    zDim=size(HessMat,3);
    FIM=zeros(xDim,xDim);
    
    C=zeros(xDim,xDim);
    for curMeas=1:numZ
        if(size(RInv,3)>1)
            RInvCur=RInv(:,:,curMeas);
        else
            RInvCur=RInv; 
        end
        FIM=FIM+JacobMat'*RInv*JacobMat;
        
        RhzVal=RInvCur*(h-z(:,curMeas));
        for k=1:xDim
            C(:,k)=reshape(HessMat(:,k,:),[xDim,zDim])*RhzVal;
        end
        FIM=FIM+C;
    end
else
    numMeas=size(RInv,3);
    xDim=size(JacobMat,2);
    FIM=zeros(xDim,xDim);
    for k=1:numMeas
        FIM=FIM+JacobMat'*RInv(:,:,k)*JacobMat;
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
