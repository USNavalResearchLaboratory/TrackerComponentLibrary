function [xIntervalEst,PIntervalEst,xFwdPred,PFwdPred,xFwdPost,PFwdPost]=KalmanIntervalSmoother(xFwdPred,PFwdPred,xFwdPost,PFwdPost,N,zCur,RCur,HCur,FInterval,QPrev,uPrev,hasLastPred,hasLastUpdate)
%%KALMANINTERVALSMOOTHER Update the standard forward-backward Kalman
%                smoother over a sliding or growing interval given a new
%                measurement. The smoothed estimates over the entire
%                interval are available. If one already computed the
%                predicted state at the current time and the updated state
%                at the current time, they can be provided rather than
%                being recomputed. This can slide an interval of estimates
%                fixed at length N or grow an interval until it is length
%                N. See the example below for how one can grow the interval
%                until it is a desired length.
%
%INPUTS: xFwdPred The xDimX(N-1) set of predicted target states. If the
%                 interval is is increasing in size (has not yet reached
%                 the desired N), then the last dimension is NCur-2 in
%                 size, where NCur is the current length, unless
%                 hasLastPred=true, in which case it is NCur-1 in length
%                 with the final prediction already computed. If 
%                 hasLastPred=false, then the final column of xFwdPred is
%                 the state predicted to the time-step PRIOR to the current
%                 one. Otherwise, if hasLastPred=true, then the final
%                 column is the state predicted to the time of the current
%                 measurement zCur. 
%        PFwdPred The xDimXxDimX(N-1) set of predicted state covariance
%                 matrices associated with xFwdPred. The size and meaning
%                 of the final dimension is the same as for xFwdPred.
%        xFwdPost The xDimXN set of posterior state estimates. If the
%                 interval is is increasing in size (has not yet reached
%                 the desired N), then the last dimension is NCur-1 in
%                 size, where NCur is the current length, unless
%                 hasLastUpdate=true, in which case it is NCur in length
%                 with the final posterior state already computed.  If 
%                 hasLastUpdate=false, then the final column of xFwdPost is
%                 the posterior state at the time-step PRIOR to the current
%                 one. Otherwise, if hasLastUpdate=true, then the final
%                 column is the posterior state at the time of the current
%                 measurement zCur. 
%        PFwdPost The xDimXxDimX(N-1) set of posterior state covariance
%                 matrices associated with xFwdPost. The size and meaning
%                 of the final dimension is the same as for xFwdPost.
%               N The positive scalar integer desired length of the
%                 smoothed interval; N>=2. If omitted or an empty matrix is
%                 passed, then N is set to the number of columns of
%                 xFwdPost.
%            zCur The zDimX1 measurement at the current time. If
%                 hasLastUpdate=true, then an empty matrix can be passed.
%            RCur The zDimXzDim measurement covariance matrix associated
%                 with zCur. If hasLastUpdate=true, then an empty matrix
%                 can be passed.
%            HCur The zDimXxDim measurement matrix. If hasLastUpdate=true,
%                 then an empty matrix can be passed.
%       FInterval The xDimXxDimX(N-1) set of xDimXxDim state prediction
%                 matrices over the entire region. If the interval is still
%                 increasing in length(N>NCur), then the last dimension is
%                 NCur. FInterval(:,:,end) predicts the posterior of the
%                 previous time-step to the current time-step. If a single
%                 xDimXxDim matrix is passed, then it is assumed that the
%                 state transition matrix is the same over the entire
%                 interval.
%           QPrev The xDimXxDim state prediction covariance matrix to
%                 predict the target state from the previous step to the
%                 time when zCur is obtained. If hasLastPred=true, then an
%                 empty matrix can be passed.
%           uPrev The xDimX1 control input to the state prediction. If
%                 there is no input or hasLastPred=true, then an empty
%                 matrix can be passed or this input can be omitted.
%     hasLastPred This boolean value indicates whether xFwdPred and
%                 PFwdPred contain the prediction to the time of zCur. The
%                 default if omitted or an empty matrix is passed is false.
%   hasLastUpdate This boolean value indicates whether xFwdPost and
%                 xFwdPost contain the posterior values at the time of
%                 zCur. The default if omitted or an empty matrix is passed
%                 is false.
%
%OUTPUTS: xIntervalEst The xDimXN set of smoothed state estimates over the
%                      interval. The oldest estimate is xIntervalEst(:,1)
%                      and the newest xIntervalEst(:,N). If the interval is
%                      expanding, the second dimension will be size NCur.
%         PIntervalEst The xDimXxDimXN set of covariance matrices
%                      associated with those in xIntervalEst. If the
%                      interval is expanding, the third dimension will be
%                      size NCur.
%    xFwdPred,PFwdPred The xDimX(N-1) set of predicted target state and the
%                      associated xDimXxDimX(N-1) set of predicted target
%                      state covariance matrices. If the interval is
%                      expanding, the final dimensions will be length
%                      NCur-1.
%    xFwdPost,PFwdPost The xDimXN set of posterior state estimates and the
%                      associated xDimXxDimXN set of covariance matrices.
%                      If the interval is expanding, the final dimensions
%                      will be length NCur.
% 
%This function implements the interval-based Kalman smoother that is
%described in Chapter 8.6 of [1].
%
%EXAMPLE:
%We start a track with one-point initialization. The results are then fed
%to a Kalman interval smoother until the interval reaches length N. After
%that, the interval slides. We also consider a Kalman filter and show that
%it produces identical estimates to those at the front of the interval
%smoother. The measurement error reduction factor (MERF) and the NEES of
%the estimates of the Kalman filter and the Kalman smoother are compared.
%The MERF of the smoother is lower (better) and the NEES of both are close
%to 1, indicating estimator consistency.
% T=1;%Sample period.
% numRuns=500;
% N=4;%The length of the interval considered.
% numSteps=100;
% 
% %Parameters for the noise process dynamics and measurement.
% H=[eye(2,2),zeros(2,2)];
% zDim=size(H,1);
% xDim=size(H,2);
% R=diag([40;40]);
% SR=chol(R,'lower');
% 
% %Statistics for the initial state.
% x0Mean=[1e3;0;75;-50];
% x0Cov=diag([1e3^2;100^2;25^2;25^2]);
% x0S=chol(x0Cov,'lower');
% 
% %Parameters for the state dynamics.
% q=processNoiseSuggest('PolyKal-ROT',9.8,1);
% F=FPolyKal(T,xDim,1);
% Q=QPolyKal(T,xDim,1,q);
% SQ=chol(Q,'lower');
% maxVel=400;%Used for single-point initiation.
% higherDerivStdDev=maxVel/sqrt(2);
% 
% xTrue=zeros(xDim,numSteps+N-1,numRuns);
% z=zeros(zDim,numSteps+N-1,numRuns);
% xSmoothed=zeros(xDim,numSteps,numRuns);
% PSmoothed=zeros(xDim,xDim,numSteps,numRuns);
% xKalman=zeros(xDim,numSteps,numRuns);
% PKalman=zeros(xDim,xDim,numSteps,numRuns);
% xKalmanAlt=zeros(xDim,numSteps,numRuns);
% PKalmanAlt=zeros(xDim,xDim,numSteps,numRuns);
% for curRun=1:numRuns
%     %Draw the initial state.
%     xTrue(:,1,curRun)=x0Mean+x0S*randn(xDim,1);
%     %Get the first measurement.
%     z(:,1,curRun)=H*xTrue(:,1,curRun)+SR*randn(zDim,1);
% 
%     %Get the states and measurements for the rest of the steps plus extra
%     %steps so that the Kalman filter and smoother cover the same range.
%     for curStep=2:(numSteps+N-1)
%         xTrue(:,curStep,curRun)=F*xTrue(:,curStep-1,curRun)+SQ*randn(xDim,1);
%         z(:,curStep,curRun)=H*xTrue(:,curStep,curRun)+SR*randn(zDim,1);
%     end
%     
%     %Start the filters with one-point initiation.
%     [xInit,PInit]=onePointCartInit(z(:,1,curRun),SR,higherDerivStdDev);
%     
%     %Initialize the Kalman smoother.
%     xFwdPred=[];
%     PFwdPred=[];
%     xFwdPost=xInit;
%     PFwdPost=PInit;
%     
%     %Initialize the Kalman filter.
%     xKalman(:,1,curRun)=xInit;
%     PKalman(:,:,1,curRun)=PInit;
%     xKalmanAlt(:,1,curRun)=xInit;
%     PKalmanAlt(:,:,1,curRun)=PInit;
%     for curStep=2:numSteps
%         [xIntervalEst,PIntervalEst,xFwdPred,PFwdPred,xFwdPost,PFwdPost]=KalmanIntervalSmoother(xFwdPred,PFwdPred,xFwdPost,PFwdPost,N,z(:,curStep,curRun),R,H,F,Q);
%         if(curStep>=N)
%             xSmoothed(:,curStep-N+1,curRun)=xIntervalEst(:,1);
%             PSmoothed(:,:,curStep-N+1,curRun)=PIntervalEst(:,:,1);
%         end
%         xKalmanAlt(:,curStep,curRun)=xIntervalEst(:,end);
%         PKalmanAlt(:,:,curStep,curRun)=PIntervalEst(:,:,end);
%         
%         [xPred,PPred]=discKalPred(xKalman(:,curStep-1,curRun),PKalman(:,:,curStep-1,curRun),F,Q);
%         [xUpdate,PUpdate]=KalmanUpdate(xPred,PPred,z(:,curStep,curRun),R,H);
%         xKalman(:,curStep,curRun)=xUpdate;
%         PKalman(:,:,curStep,curRun)=PUpdate;
%     end
%     
%     %Now, go N-1 steps more with the smoother, so that the Kalman filter
%     %estimates and the smoother estimates cover the same range of values.
%     for curStep=(numSteps+1):(numSteps+N-1)
%         [xIntervalEst,PIntervalEst,xFwdPred,PFwdPred,xFwdPost,PFwdPost]=KalmanIntervalSmoother(xFwdPred,PFwdPred,xFwdPost,PFwdPost,N,z(:,curStep,curRun),R,H,F,Q);
%         xSmoothed(:,curStep-N+1,curRun)=xIntervalEst(:,1);
%         PSmoothed(:,:,curStep-N+1,curRun)=PIntervalEst(:,:,1);
%     end
% end
% 
% %The maximum absolute difference between the standard Kalman filter
% %estimate and the estimate at the front of the Kalman smoother interval.
% %This should be 0.
% maxKalDiff=max(abs(xKalman(:)-xKalmanAlt(:)))
% 
% xTrue=xTrue(:,1:numSteps,:);
% z=z(:,1:numSteps,:);
% posTrue=xTrue(1:2,1:numSteps,:);
% 
% MERFSmoothed=calcMERF(posTrue,z,xSmoothed(1:2,:,:),0,true);
% MERFKalman=calcMERF(posTrue,z,xKalman(1:2,:,:),0,true);
% NEESSmoothed=calcNEES(xTrue,xSmoothed,PSmoothed);
% NEESKalman=calcNEES(xTrue,xKalman,PKalman);
% 
% figure(1)
% clf
% hold on
% plot(MERFKalman,'-k','linewidth',4)
% plot(MERFSmoothed,'--r','linewidth',2)
% h1=xlabel('Step');
% h2=ylabel('MERF');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Kalman Filter','Kalman Smoother')
% 
% figure(2)
% clf
% hold on
% plot(NEESKalman,'-k','linewidth',4)
% plot(NEESSmoothed,'--r','linewidth',2)
% h1=xlabel('Step');
% h2=ylabel('NEES');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Kalman Filter','Kalman Smoother','location','southeast')
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(HCur,2);
NCur=size(xFwdPost,2);

if(nargin<13||isempty(hasLastUpdate))
    hasLastUpdate=false;
end

if(nargin<12||isempty(hasLastPred))
    hasLastPred=false;
end

if(nargin<11)
    uPrev=[];
end

if(isempty(N))
    %Assume that the batch is already the desired size if the batch length
    %is not specified.
    N=NCur;
end

if(N>NCur)
    NFEnd=NCur;
else
    NFEnd=N-1;
end

if(size(FInterval,3)==1)
    FInterval=repmat(FInterval,[1,1,NFEnd]);
end

if(hasLastPred==false)    
    xFwdPostEnd=xFwdPost(:,NCur);
    PFwdPostEnd=PFwdPost(:,:,NCur);
    [xFwdPredEnd,PFwdPredEnd]=discKalPred(xFwdPostEnd,PFwdPostEnd,FInterval(:,:,NFEnd),QPrev,uPrev);

    if(NCur==N)
        %If the interval is sliding.
        xFwdPred=[xFwdPred(:,2:(N-1)),xFwdPredEnd];
        PFwdPred=cat(3,PFwdPred(:,:,2:(N-1)),PFwdPredEnd);
    else%The interval is growing. NCur should be  the length of xFwdPred+1
        %if the interval is growing.
        xFwdPred(:,NCur)=xFwdPredEnd;
        PFwdPred(:,:,NCur)=PFwdPredEnd;
    end
    hasLastUpdate=false;
end

if(hasLastUpdate==false)
    [xFwdEnd,PFwdEnd]=KalmanUpdate(xFwdPred(:,NFEnd),PFwdPred(:,:,NFEnd),zCur,RCur,HCur);

    if(NCur==N)
        %If the interval is sliding.
        xFwdPost=[xFwdPost(:,2:N),xFwdEnd];
        PFwdPost=cat(3,PFwdPost(:,:,2:N),PFwdEnd);
    else%The interval is growing.
        xFwdPost(:,end+1)=xFwdEnd;
        PFwdPost(:,:,end+1)=PFwdEnd;
        NCur=NCur+1;
    end
end

%Run the backwards Kalman smoother
xIntervalEst=zeros(xDim,NCur);
PIntervalEst=zeros(xDim,xDim,NCur);
xIntervalEst(:,NCur)=xFwdPost(:,NCur);
PIntervalEst(:,:,NCur)=PFwdPost(:,:,NCur);
for curStep=(NCur-1):-1:1
    C=PFwdPost(:,:,curStep)*FInterval(:,:,curStep)'/PFwdPred(:,:,curStep);
    xIntervalEst(:,curStep)=xFwdPost(:,curStep)+C*(xIntervalEst(:,curStep+1)-xFwdPred(:,curStep));
    PIntervalEst(:,:,curStep)=PFwdPost(:,:,curStep)+C*(PIntervalEst(:,:,curStep+1)-PFwdPred(:,:,curStep))*C';
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
