function [yIntervalEst,PInvIntervalEst,yFwdPred,PInvFwdPred,yFwdEnd,PInvFwdEnd]=FPInfoIntervalSmoother(yFwdPred,PInvFwdPred,yFwdPrev,PInvFwdPrev,N,z,R,H,F,Q,u,hasLastPred)
%%FPINFOINTERVALSMOOTHER Update the Fraser-Potter information smoother 
%                smoother over a sliding or growing interval given a new
%                measurement. The smoothed estimates over the entire
%                interval are available. If one already computed the
%                predicted state at the current time, it can be provided
%                rather than being recomputed. This can slide an interval
%                of estimates fixed at length N or grow an interval until
%                it is length N. See the example below for how one can grow
%                the interval until it is a desired length.
%
%INPUTS: yFwdPred The xDimXN set of predicted information states. If the
%                 interval is is increasing in size (has not yet reached
%                 the desired N), then the last dimension is NCur-1 in
%                 size, where NCur is the current length of z, unless
%                 hasLastPred=true, in which case it is NCur in length
%                 with the final prediction already computed. If 
%                 hasLastPred=false, then the final column of xFwdPred is
%                 the state predicted to the time-step PRIOR to the current
%                 one. Otherwise, if hasLastPred=true, then the final
%                 column is the state predicted to the time of the latest
%                 measurement in z. This cannot be an empty matrix.
%     PInvFwdPred The xDimXxDimXN set of predicted inverse covariance
%                 matrices associated with yFwdPred. The size and meaning
%                 of the final dimension is the same as for yFwdPred.
%        yFwdPrev The xDimX1 posterior information state estimate at the
%                 time of the second to latest measurement in z. When
%                 calling this function repeatedly, this can be the yFwdEnd
%                 output from the last function call. If hasLastPred=true,
%                 then an empty matrix can be passed in place of this
%                 input.
%     PInvFwdPrev The xDimXxDim inverse covariance matrix associate with
%                 yFwdPrev. When calling this function repeatedly, this can
%                 be the PInvFwdEnd output from the last function call. If
%                 hasLastPred=true, then an empty matrix can be passed in
%                 place of this input.
%               N The positive scalar integer desired length of the
%                 smoothed interval; N>=2. If omitted or an empty matrix is
%                 passed, then N is set to the number of columns of z
%                 (NCur).
%               z The zDimXN set of measurements over the interval. z(:,N)
%                 is the latest measurement. If the interval is expanding,
%                 then z is zDimXNCur in length.
%               R The zDimXzDimXN set of measurement covariance matrices
%                 associated with the values in z. If the interval is
%                 expanding, then the last dimension is NCur in size. If
%                 all the matrices are the same, then a single zDimXzDim
%                 matrix can be passed.
%               H The xDimXxDimXN set of measurement matrices associated
%                 with the values in z. If the interval is expanding, then
%                 the last dimension is NCur in size. If all the matrices
%                 are the same, then a single xDimXzDim matrix can be
%                 passed.
%               F The xDimXxDimX(N-1) set of state transition matrices over
%                 the interval. F(:,:,N-1) is the state transition matrix
%                 to the time of the latest measurerment in z. If the batch
%                 is expanding, then the last dimension is NCur-1 in size.
%                 If all F are the same, then a single xDimXxDim matrix can
%                 be passed.
%               Q The xDimXxDimX(N-1) set of process noise covariance
%                 matrices. the meaning of the third dimension is the same
%                 as for F. If all Q are the same, then a single xDimXxDim
%                 matrix can be passed.
%               u The xDimX(N-1) set of state transition control inputs,
%                 where the second dimension has the same meaning as the
%                 third dimensions of F. . If an empty matrix is passed,
%                 then it is assumed that all control inputs are zero. If
%                 aall the u vectors are the same, then a single xDimX1
%                 vector can be passed.
%     hasLastPred This boolean value indicates whether yFwdPred and
%                 PInvFwdPred contain the prediction to the time of the
%                 latest measurement in z. The default if omitted or an
%                 empty matrix is passed is false.
%
%OUTPUTS: yIntervalEst The xDimXN set of smoothed information state
%                      estimates over the interval. The oldest estimate is
%                      yIntervalEst(:,1)
%                      and the newest yIntervalEst(:,N). If the interval is
%                      expanding, the second dimension will be size NCur.
%      PInvIntervalEst The xDimXxDimXN set of inverse covariance matrices
%                      associated with those in yIntervalEst. If the
%                      interval is expanding, the third dimension will be
%                      size NCur.
% yFwdPred,PInvFwdPred The xDimXN set of predicted information states and
%                      the associated xDimXxDimX(N-1) set of predicted
%                      inverse covariance matrices. If the interval is
%                      expanding, the final dimensions will be length
%                      NCur.
%   yFwdEnd,PInvFwdEnd The xDimX1 and xDimXxDim posterior information state
%                      and its inverse covariance matrix at the latest
%                      time.
%
%The Kalman smoothing technique utilizing forward and backwards filters is
%described in [1]. The use with information filters is related.
%
%EXAMPLE:
%Given an uninformative prior, the information smoother is initialized with
%the output of an information filter at the first step. Then, as more
%measurements are obtained, the batch length increases until it reached
%length N. After that, the interval slides. We also run an information
%filter and show that not only does the smoothed estimate at the end of the
%interval (the most smoothed) outperform the information filter, but that
%the estimate at the start of the interval is the same as the real-time
%information filter output. The measurement error reduction factor (MERF)
%and the NEES of the estimates (converted to target states) of the
%information filter and the information  smoother are compared. The MERF of
%the smoother is lower (better) and the NEES of both are close to 1,
%indicating estimator consistency.
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
% 
% %The initial uninformative state.
% yInit=zeros(xDim,1);
% PInvInit=zeros(xDim,xDim,1);
% 
% xTrue=zeros(xDim,numSteps+N-1,numRuns);
% z=zeros(zDim,numSteps+N-1,numRuns);
% ySmoothed=zeros(xDim,numSteps,numRuns);
% PInvSmoothed=zeros(xDim,xDim,numSteps,numRuns);
% yInfo=zeros(xDim,numSteps,numRuns);
% PInvInfo=zeros(xDim,xDim,numSteps,numRuns);
% yInfoAlt=zeros(xDim,numSteps,numRuns);
% PInvInfoAlt=zeros(xDim,xDim,numSteps,numRuns);
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
%     %Initialize the information filters and the information smoother with a
%     %single measurement.
%     [yUpdate,PInvUpdate]=infoFilterUpdate(yInit,PInvInit,z(:,1,curRun),R,H);
%     yFwdPred=yInit;
%     PInvFwdPred=PInvInit;
%     yFwdPrev=yUpdate;
%     PInvFwdPrev=PInvUpdate;
%     
%     yInfo(:,1,curRun)=yUpdate;
%     PInvInfo(:,:,1,curRun)=PInvUpdate;
%     
%     yInfoAlt(:,1,curRun)=yUpdate;
%     PInvInfoAlt(:,:,1,curRun)=PInvUpdate;
%     
%     %Continue with the rest of the measurements.
%     for curStep=2:numSteps
%         zSpan=z(:,max(curStep-N+1,1):curStep,curRun);
%         
%         [yIntervalEst,PInvIntervalEst,yFwdPred,PInvFwdPred,yFwdPrev,PInvFwdPrev]=FPInfoIntervalSmoother(yFwdPred,PInvFwdPred,yFwdPrev,PInvFwdPrev,N,zSpan,R,H,F,Q);
%         if(curStep>=N)
%             ySmoothed(:,curStep-N+1,curRun)=yIntervalEst(:,1);
%             PInvSmoothed(:,:,curStep-N+1,curRun)=PInvIntervalEst(:,:,1);
%         end
%         yInfoAlt(:,curStep,curRun)=yIntervalEst(:,end);
%         PInvInfoAlt(:,:,curStep,curRun)=PInvIntervalEst(:,:,end);
%         
%         [yPred,PInvPred]=infoFilterDiscPred(yInfo(:,curStep-1,curRun),PInvInfo(:,:,curStep-1,curRun),F,Q);
%         [yUpdate,PInvUpdate]=infoFilterUpdate(yPred,PInvPred,z(:,curStep,curRun),R,H);
%         yInfo(:,curStep,curRun)=yUpdate;
%         PInvInfo(:,:,curStep,curRun)=PInvUpdate;
%     end
%     
%     %Now, go N-1 steps more with the smoother, so that the information
%     %filter estimates and the information smoother estimates cover the same
%     %range of values.
%     for curStep=(numSteps+1):(numSteps+N-1)
%         zSpan=z(:,(curStep-N+1):curStep,curRun);
%         [yIntervalEst,PInvIntervalEst,yFwdPred,PInvFwdPred,yFwdPrev,PInvFwdPrev]=FPInfoIntervalSmoother(yFwdPred,PInvFwdPred,yFwdPrev,PInvFwdPrev,N,zSpan,R,H,F,Q);
%         ySmoothed(:,curStep-N+1,curRun)=yIntervalEst(:,1);
%         PInvSmoothed(:,:,curStep-N+1,curRun)=PInvIntervalEst(:,:,1);
%     end
% end
% 
% %The maximum absolute difference between the standard information filter
% %estimate and the estimate at the front of the infomration smoother batch.
% %This should be 0.
% maxInfoDiff=max(abs(yInfo(:)-yInfoAlt(:)))
% 
% %Get the target states. We shall skip the first one which is not completely
% %observable.
% xSmoothed=zeros(xDim,numSteps-1,numRuns);
% PSmoothed=zeros(xDim,xDim,numSteps-1,numRuns);
% xEst=zeros(xDim,numSteps-1,numRuns);
% PEst=zeros(xDim,xDim,numSteps-1,numRuns);
% for curRun=1:numRuns 
%     for curStep=2:numSteps
%         xSmoothed(:,curStep-1,curRun)=PInvSmoothed(:,:,curStep,curRun)\ySmoothed(:,curStep,curRun);
%         PSmoothed(:,:,curStep-1,curRun)=inv(PInvSmoothed(:,:,curStep,curRun));
%         
%         xEst(:,curStep-1,curRun)=PInvInfo(:,:,curStep,curRun)\yInfo(:,curStep,curRun);
%         PEst(:,:,curStep-1,curRun)=inv(PInvInfo(:,:,curStep,curRun));
%     end
% end
% 
% %Limit the period plotted to the range where the target state is fully
% %observable and is common to the smoother and the filter.
% xTrue=xTrue(:,2:numSteps,:);
% z=z(:,2:numSteps,:);
% posTrue=xTrue(1:2,:,:);
% 
% MERFSmoothed=calcMERF(posTrue,z,xSmoothed(1:2,:,:),0,true);
% MERFInfo=calcMERF(posTrue,z,xEst(1:2,:,:),0,true);
% NEESSmoothed=calcNEES(xTrue,xSmoothed,PSmoothed);
% NEESInfo=calcNEES(xTrue,xEst,PEst);
% 
% figure(1)
% clf
% hold on
% plot(2:numSteps,MERFInfo,'-k','linewidth',4)
% plot(2:numSteps,MERFSmoothed,'--r','linewidth',2)
% h1=xlabel('Step');
% h2=ylabel('MERF');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Information Filter','Information Smoother')
% 
% figure(2)
% clf
% hold on
% plot(2:numSteps,NEESInfo,'-k','linewidth',4)
% plot(2:numSteps,NEESSmoothed,'--r','linewidth',2)
% h1=xlabel('Step');
% h2=ylabel('NEES');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% legend('Information Filter','Information Smoother','location','southeast')
%
%REFERENCES:
%[1] D. C. Fraser and J. E. Potter, "The optimal linear smoother as a 
%    combination of two optimum linear filters," IEEE Transactions on
%    Automatic Control, pp. 387-390, Aug. 1969.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(yFwdPred,1);
NCur=size(z,2);

if(nargin<12||isempty(hasLastPred))
    hasLastPred=false;
end

if(isempty(N))
    %Assume that the batch is already the desired size if the batch length
    %is not specified.
    N=NCur;
end

NFEnd=NCur-1;

if(size(F,3)==1)
    F=repmat(F,[1,1,NFEnd]);
end

if(size(Q,3)==1)
    Q=repmat(Q,[1,1,NFEnd]);
end

if(nargin<11||isempty(u))
    u=zeros(xDim,NFEnd);
elseif(size(u,2)==1)
    u=repmat(u,[1,NFEnd]);
end

if(size(H,3)==1)
    H=repmat(H,[1,1,NCur]); 
end

if(size(R,3)==1)
    R=repmat(R,[1,1,NCur]); 
end

if(hasLastPred==false)
    [yFwdPredEnd,PInvFwdPredEnd]=infoFilterDiscPred(yFwdPrev,PInvFwdPrev,F(:,:,NFEnd),Q(:,:,NFEnd),u(:,NFEnd));

    if(NCur==N&&size(yFwdPred,2)==N)
        %If the interval is sliding.
        yFwdPred=[yFwdPred(:,2:N),yFwdPredEnd];
        PInvFwdPred=cat(3,PInvFwdPred(:,:,2:N),PInvFwdPredEnd);
    else%The interval is growing; NCur is one more than the current length
        %of yFwdPred.
        yFwdPred(:,NCur)=yFwdPredEnd;
        PInvFwdPred(:,:,NCur)=PInvFwdPredEnd;
    end
end

%Run the information filter backwards.
yRev=zeros(xDim,NCur);
PInvRev=zeros(xDim,xDim,NCur);

[yRev(:,NCur),PInvRev(:,:,NCur)]=infoFilterUpdate(zeros(xDim,1),zeros(xDim,xDim),z(:,NCur),R(:,:,NCur),H(:,:,NCur));
for curStep=(NCur-1):-1:1
    [yRevPred,PInvRevPred]=infoFilterDiscPredRev(yRev(:,curStep+1),PInvRev(:,:,curStep+1),F(:,:,curStep),Q(:,:,curStep),u(:,curStep));
    [yRev(:,curStep),PInvRev(:,:,curStep)]=infoFilterUpdate(yRevPred,PInvRevPred,z(:,curStep),R(:,:,curStep),H(:,:,curStep));
end

PInvIntervalEst=PInvFwdPred+PInvRev;
yIntervalEst=yFwdPred+yRev;

if(nargout>4)
    yFwdEnd=yIntervalEst(:,end);
    PInvFwdEnd=PInvIntervalEst(:,:,end);
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
