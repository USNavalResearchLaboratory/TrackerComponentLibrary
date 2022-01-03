function [xUpdate,PUpdate,innov,Pzz,W]=KalmanUpdate(xPred,PPred,z,R,H,condIssue)
%KALMANUPDATE Perform the measurement update step in the standard linear
%             Kalman filter. This assumes that the predicted target state
%             and the measurement are both multivariate Gaussian
%             distributed.
%
%INPUTS: xPred The xDimX1 predicted target state.
%        PPred The xDimXxDim predicted state covariance matrix.
%            z The zDimX1 measurement vector.
%            R The zDimXzDim measurement covariance matrix. 
%            H The optional zDimXxDim measurement matrix. The measurement
%              is modeled as z=H*x+noise. If this parameter is omitted or
%              an empty matrix is passed, then H will be taken as a
%              zDimXzDim identity matrix followed by columns of zeros
%              (Assuming that zDim<=xDim. Otherwise, H must be provided).
%    condIssue An optional parameter indicating whether a conditioning
%              issue with the measurement covariance matrix. Singularity
%              can occur if two elements are perfectly correlated. This
%              replaces the matrix inverse with a pseudoinverse.
%
%OUTPUTS: xUpdate The xDim X 1 updated target state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%Given a prediction of the target state xPred with an associated covariance
%matrix PPred, and assuming that they are the mean and covariance matrix of
%a Gaussian distribution, under the measurement model
%z=H*x+w
%where z is the measurement, x is the true target state, H is a zDimXxDim
%matrix and w is additive Gaussian noise, this function updates the state
%estimate and its associated covariance estimate.
%
%The Joseph-form covariance update is used for improved numerical
%stability. The algorithm is derived in Chapter 5 of [1].
%
%EXAMPLE:
%This example simulates a linear dynamic process that extracts Cartesian
%measurements. it then looks at the measurement error reduction factor
%(MERF), and the normalized estimation error squared (NEES). The MERF
%should be less than 1 or else there is no point in using the filter-
%connecting the dots works better. The NEES should be close to 1 to
%indicate covariance consistency.
% numMonteCarlo=500;
% numSteps=50;
% %Power spectral density of the process noise.
% q0=processNoiseSuggest('PolyKal-ROT',2*9.8,1);
% xDim=4;%2D position and velocity.
% T=1;%prediction interval.
% F=FPolyKal(T,xDim,1);
% Q=QPolyKal(T,xDim,1,q0);
% SQ=chol(Q,'lower');
% R=[50^2,-1000;
%    -1000,50^2];
% SR=chol(R,'lower');
% H=[eye(2,2),zeros(2,2)];
% zDim=size(H,1);
% 
% z=zeros(zDim,numSteps,numMonteCarlo);
% xTrue=zeros(xDim,numSteps,numMonteCarlo);
% %The initial state.
% xTrue(:,1,:)=repmat([1e3;0;75;-50],[1,1,numMonteCarlo]);
% maxVel=400;%Used for single-point initiation.
% higherDerivStdDev=maxVel/sqrt(2);%used for single-point initiation.
% 
% xEst=zeros(xDim,numSteps,numMonteCarlo);
% PEst=zeros(xDim,xDim,numSteps,numMonteCarlo);
% for curRun=1:numMonteCarlo
%     z(:,1,curRun)=H*xTrue(:,1,curRun)+SR*randn(zDim,1);
%     for curStep=2:numSteps        
%         xTrue(:,curStep,curRun)=F*xTrue(:,curStep-1,curRun)+SQ*randn(xDim,1);
%         z(:,curStep,curRun)=H*xTrue(:,curStep,curRun)+SR*randn(zDim,1);
%     end
% 
%     %One-point initialization. This works in this instance, because H
%     %extracts the Cartesian components.
%     [xInit,PInit]=onePointCartInit(z(:,1,curRun),R,higherDerivStdDev,1);
%     xEst(:,1,curRun)=xInit;
%     PEst(:,:,1,curRun)=PInit;
%     for curStep=2:numSteps
%         [xPred,PPred]=discKalPred(xEst(:,curStep-1,curRun,1),PEst(:,:,curStep-1,curRun,1),F,Q);
%         [xUpdate,PUpdate]=KalmanUpdate(xPred,PPred,z(:,curStep,curRun),R,H);
%         xEst(:,curStep,curRun)=xUpdate;
%         PEst(:,:,curStep,curRun)=PUpdate;
%     end
% end
% 
% posTrue=xTrue(1:2,:,:);
% %The MERF requires the measurement to be in the same coordinates are the
% %corresponding element sof the state.
% MERF=calcMERF(posTrue,z,xEst(1:2,:,:),0,true);
% NEES=calcNEES(xTrue,xEst,PEst);
% figure(1)
% clf
% hold on
% plot(MERF,'-k','linewidth',2)
% h1=xlabel('step');
% h2=ylabel('Measurement Error Reduction Factor');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% figure(2)
% clf
% hold on
% plot(NEES,'-r','linewidth',2)
% h1=xlabel('step');
% h2=ylabel('NEES');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPred,1);

if(nargin<6||isempty(condIssue))
    condIssue=false;
end

if(nargin<5||isempty(H))
	zDim=size(z,1); 
    H=[eye(zDim,zDim),zeros(zDim,xDim-zDim)];
end

zPred=H*xPred;
innov=z-zPred;

Pzz=R+H*PPred*H';
%Ensure symmetry
Pzz=(Pzz+Pzz')/2;

if(condIssue==false)
    opts.SYM=true;
    opts.RECT=false;
    opts.TRANSA=false;
    W=linsolve(Pzz,H*PPred,opts)';%W=PPred*H'/Pzz;
else
    %If there is possibly poor conditioning.
    W=PPred*H'*pinv(Pzz);
end

xUpdate=xPred+W*innov;

temp=W*H;
temp=eye(xDim,xDim)-temp;
PUpdate=temp*PPred*temp'+W*R*W';

%Ensure symmetry
PUpdate=(PUpdate+PUpdate')/2;
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
