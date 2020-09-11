function [xEst,PEst]=batchLSLinMeasLinDyn(z,H,F,R,kD,Q,numMeas)
%%BATCHLSLINMEASLINDYN Perform batch least squares state estimation under a
%                linear measurement model and a linear dynamic model with
%                optional process noise.
%
%INPUTS: z The zDim X N matrix of measurements for the whole batch. It is
%          assumed that the measurements have the same dimensionality over
%          the batch. If an empty matrix is passed, then it is assumed that
%          the user only wants the covariance matrix and not xEst. 
%        H The zDim X xDim X N hypermatrix of measurement matrices such
%          that H(:,:,k)*x+w is the measurement at time k, where x is the
%          state and w is zero-mean Gaussian noise with covariance matrix
%          R(:,:,k). Alternatively, if all of the measurement matrices are
%          the same, one can just pass a single zDim X xDim matrix.
%        F An xDim X xDim X (N-1) hypermatrix of  matrices. The state at
%          discrete-time k+1 is modeled as F(:,:,k) times the state at time
%          k plus zero-mean Gaussian process noise with covariance matrix
%          Q(:,:,k). Alternatively, if all of the state transition matrices
%          are the same, one can just pass a single xDim X xDim matrix.
%          Note that all of the F matrices must be invertible if kD is not
%          equal to one.
%        R The zDim X zDim X N hypermatrix of measurement covariance
%          matrices. Alternatively, if all of the measurement covariance
%          matrices are the same, one can just pass a single zDimXzDim
%          matrix.
%       kD The discrete time-step for which the covariance matrix of an
%          estimate is desired, where z(:,1) is at discrete time-step 1
%          (not 0).
%        Q The xDim X xDim X (N-1) hypermatrix of process noise covariance
%          matrices. Alternatively, if all of the process noise covariance
%          matrices are the same, one can just pass a single xDimXxDim
%          matrix. If an empty matrix is passed for Q, then the estimation
%          is performed assuming there is no process noise.
%  numMeas If an empty matrix is passed for z, then the numMeas parameter
%          must be passed indicating the number of measurements in the
%          batch. Otherwise, this parameter is
%                ignored.
%
%OUTPUTS: xEst The batch state estimate at step kD, unless an empty matrix
%              was passed for z, in which case xEst is empty.
%         PEst A covariance matrix estimate that goes with the state
%              estimate.
%
%The algorithm is an implementation of the method of Section 3.3.2 of [1]
%for a linear dynamic model. Note that the cases in equation 28 to 32 of
%the paper did not cover all of the possibilities, so the correct
%generalization had to be derived from E[epsilon_{p,k}*epsilon_{q,k}'].
%
%REFERENCES:
%[1] A. B. Poore, B. J. Slocumb, B. J. Suchomel, F. H. Obermeyer, S. M.
%    Herman, and S. M. Gadaleta, "Batch maximum likelihood (ML) and maximum
%    a posteriori (MAP) estimation with process noise for tracking
%    applications," in Proceedings of SPIE: Signal and Data Processing of
%    Small Targets, vol. 5204, San Diego, CA, 3 Aug. 2003, pp. 188-199.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isempty(z))
    numMeas=size(z,2);
end

xDim=size(F,1);
zDim=size(R,1);

if(size(H,3)==1)
    HMats=repmat(H,[1,1,numMeas]);
else
    HMats=H;
end

if(size(F,3)==1)
    F=repmat(F,[1,1,numMeas-1]);
end

if(size(R,3)==1)
    R=repmat(R,[1,1,numMeas]);
end

if(isempty(Q))
    Q=zeros(xDim,xDim,numMeas);
elseif(size(Q,3)==1)
    Q=repmat(Q,[1,1,numMeas-1]);
end

transMats=getTransMats(F);

WMat=zeros(zDim*numMeas,zDim*numMeas);
for p=1:numMeas
    pIdxMin=(p-1)*zDim+1;
    pIdxMax=p*zDim;
    pSpan=pIdxMin:pIdxMax;
    for q=1:numMeas
        qIdxMin=(q-1)*zDim+1;
        qIdxMax=q*zDim;
        qSpan=qIdxMin:qIdxMax;
        if(p<kD&&q<kD)
            CumMat=zeros(xDim,xDim);
            
            for i=(max(p,q)+1):kD
                CumMat=CumMat+transMats(:,:,p,i)*Q(:,:,i-1)*transMats(:,:,q,i)';
            end
            
            WMat(pSpan,qSpan)=HMats(:,:,p)*CumMat*HMats(:,:,q)'+KDelta(p-q)*R(:,:,p);
        elseif(p>kD&&q>kD)
            CumMat=zeros(xDim,xDim);
            
            for i=(kD+1):min(p,q)
                CumMat=CumMat+transMats(:,:,p,i)*Q(:,:,i-1)*transMats(:,:,q,i)';
            end
            
            WMat(pSpan,qSpan)=KDelta(p-q)*R(:,:,p)+HMats(:,:,p)*CumMat*HMats(:,:,q)';
        elseif(p==kD&&q==kD)
            WMat(pSpan,qSpan)=R(:,:,p);
        end
    end
end

%The propagated measurement matrices.
bigH=zeros(zDim*numMeas,xDim);
for curMeas=1:numMeas
    idxMin=(curMeas-1)*zDim+1;
    idxMax=curMeas*zDim;
    
    bigH(idxMin:idxMax,:)=HMats(:,:,curMeas)*transMats(:,:,curMeas,kD);
end

WInv=inv(WMat);
PEstInv=bigH'*WInv*bigH;
PEst=inv(PEstInv);
if(isempty(z))
    xEst=[];
else
    xEst=PEstInv\bigH'*(WMat\z(:));
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

