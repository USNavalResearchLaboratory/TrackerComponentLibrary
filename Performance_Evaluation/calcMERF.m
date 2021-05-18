function val=calcMERF(zTrue,zMeas,zEst,type,is3D)
%%CALCMERF Compute the measurement error reduction factor (MERF). This is a
%          measure of how well a filter does at reducting the error of
%          measurements. If ther MERF is greater than 1, then it is
%          generally better to not use the filter at all and just "connect
%          the dots" with the measurements. This function can also be used
%          to compute an estimation error reduction factor (EERF) by
%          changing the inputs, as described below.
%
%INPUTS: zTrue The noiseless values of the measurements (or for an EERF,
%              true target state components). This is either an
%              zDimXNumSamples matrix or an zDimXNXnumSamples matrix. The
%              latter formulation is useful when the MSE over multiple
%              Monte Carlo runs of an entire length-N track is desired. In
%              the first formulation, we take N=1. Alternatively, if the
%              same true value is used for all numSamples, then zTrue can
%              just be an zDimXN matrix. Alternatively, if zTrue is the
%              same for all numSamples and all N, then just an zDimX1
%              matrix can be passed. N and numSamples are inferred from
%              zMeas.
%        zMeas An zDimXnumSamples set of measurements (or measurements
%              converted to the state domain if an EERF is desired) or a
%              zDimXNXnumSamples set of measurements (if values at N times
%              are desired).
%         zEst An zDimXnumSamples set of estimates in the measurement
%              domain (for the MERF and in the state domain for the
%              EERF) or a zDimXNXnumSamples set of measurements (if values
%              at N times are desired).
%         type An optional parameter specifying the type of MERF/EERF to
%              compute.  Possible values are
%              0 (The default if omitted or an empty matrix is passed)
%                Compute the MERF using the AEE as in Equation 14 in [1].
%              1 Compute the MERF using the RMSE as after Equation 15 in
%                [1].
%              2 Compute the MERF using a type of geometric average as
%                after Equation 16 in [1].
%         is3D An optional indicating that zEMeas is 3D. This is only used
%              is zMeas is a matrix. In such an instance, there is an
%              ambiguity whether zMeas is truly 2D, so N=1, or whether
%              numSamples=1 and zMeas is 3D with the third dimension being
%              1. If this parameter is omitted or an empty matrix is
%              passed, it is assumed that N=1.
%
%OUTPUTS: val The 1XN set of scalar values of the MERF (or EERF).
%
%The MERF and the EERF are discussed in detail in [1].
%
%EXAMPLE:
%Here, we consider the MERF when looking at the estimate of a Kalman filter
%after a number of steps when dealing with a simple Gaussian model and
%Gaussian measurements. 
% H=[1,0,0,0;
%    0,1,0,0];
% zDim=size(H,1);
% xDim=size(H,2);
% T=1;
% F=FPolyKal(T,4,1);
% q0=1e-3;
% Q=QPolyKal(T,4,1,q0);
% SQ=chol(Q,'lower');
% R=eye(2);
% SR=chol(R,'lower');
% numRuns=100;
% numSteps=20;
% xInitTrue=[1000;100;20;50];
% sigmaV=200;%Standard deviation for single-point initializtion.
% zTrue=zeros(zDim,numSteps,numRuns);
% zMeas=zeros(zDim,numSteps,numRuns);
% zEst=zeros(zDim,numSteps,numRuns);
% for curRun=1:numRuns
%     %Create the true track
%     xTrue=xInitTrue;
%     zTrue(:,1,curRun)=H*xTrue;
%     zMeas(:,1,curRun)=H*xTrue+SR*randn(zDim,1);
%     [xEst,PEst]=onePointCartInit(zMeas(:,1,curRun),SR,sigmaV);
%     zEst(:,1,curRun)=H*xEst;
%     SEst=chol(PEst,'lower');
%     for curStep=1:numSteps
%         [xEst,SEst]=sqrtDiscKalPred(xEst,SEst,F,SQ);
%         
%         xTrue=F*xTrue+SQ*randn(xDim,1);
%         zTrue(:,curStep,curRun)=H*xTrue;
%         zMeas(:,curStep,curRun)=zTrue(:,curStep,curRun)+SR*randn(zDim,1);
%         [xEst,SEst]=sqrtKalmanUpdate(xEst,SEst,zMeas(:,curStep,curRun),SR,H);
%         zEst(:,curStep,curRun)=H*xEst;
%     end
% end
% calcMERF(zTrue,zMeas,zEst,0)
% calcMERF(zTrue,zMeas,zEst,1)
% calcMERF(zTrue,zMeas,zEst,2)
% %One will see that the error is reduced usually to less than 0.5 for all
% %estimators. As a sanity check, one will see that for all estimators, if
% %the measurement is used as the estimate, then the estimators return 1
% calcMERF(zTrue,zMeas,zMeas,0)
% calcMERF(zTrue,zMeas,zMeas,1)
% calcMERF(zTrue,zMeas,zMeas,2)
%
%REFERENCES:
%[1] X. R. Li and Z. Zhao, "Measures of performance for evaluation of
%    estimators and filters," in Proceedings of SPIE: Conference on Signal
%    and Data processing of Small Targets, vol. 4473, San Diego, CA, 29
%    Jul. 2001, pp. 530-541.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(is3D))
    is3D=false; 
end

if(nargin<4||isempty(type))
    type=0; 
end

zDim=size(zMeas,1);
if(ismatrix(zMeas)&&is3D==false)
    N=1;
    numSamples=size(zMeas,2);
    zMeas=reshape(zMeas,[zDim,1,numSamples]);
    
    if(size(zTrue,2)==1)
        %If the true values are the same for all samples.
        zTrue=repmat(zTrue,[1,1,numSamples]);
    else
        zTrue=reshape(zTrue,[zDim,1,numSamples]);
    end
else
    N=size(zMeas,2);
    numSamples=size(zMeas,3);
    
    if(ismatrix(zTrue))
        if(size(zTrue,2)==1)
            %If the true values are the same for all samples and for all N.
            zTrue=repmat(zTrue,[1,N,numSamples]);
        else        
            %If the true values are the same for all samples.
            zTrue=repmat(zTrue,[1,1,numSamples]);
        end
    end
end

if(ismatrix(zEst)&&is3D==false)
    zEst=reshape(zEst,[zDim,1,numSamples]);
    zEst=repmat(zEst,[1,N,1]);
end

switch(type)
    case 0%AEE
        val=calcAEE(zTrue,zEst,is3D)./calcAEE(zTrue,zMeas,is3D);
    case 1%RMSE
        val=calcRMSE(zTrue,zEst,is3D)./calcRMSE(zTrue,zMeas,is3D);
    case 2%A type of GAE
        val=zeros(1,N);
        for k=1:N
            val(k)=exp((1/(2*numSamples))*sum(log(sum((zTrue(:,k,:)-zEst(:,k,:)).^2,1)./sum((zTrue(:,k,:)-zMeas(:,k,:)).^2,1))));
        end
    otherwise
        error('Unknown type specified')
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
