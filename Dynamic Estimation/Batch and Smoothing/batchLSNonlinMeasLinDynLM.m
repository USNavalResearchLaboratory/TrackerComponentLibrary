function [xEst,PEst,xBatchEst,exitCode]=batchLSNonlinMeasLinDynLM(xInit,z,h,F,R,kD,Q,HJacob,optimParam)
%%BATCHLSNONLINMEASLINDYNLM  Given a single measurement or a batch of
%           measurements, perform batch least squares state estimation
%           refinement under a nonlinear measurement model and a linear
%           dynamic model with optional process noise. Nonlinear least
%           squares optimization is performed using the Levenburg-Marquart
%           algorithm with a numerically estimated Jacobian.
%
%INPUTS: xInit The xDimX1 initial estimate of the target state at the time
%              of the kDth measurement or, if process noise is present,
%              one can optionally provide an xDimXN matrix of initial
%              state estimate vectors over an entire batch of measurements. 
%            z The zDim X N matrix of measurements for the whole batch. It
%              is assumed that the measurements have the same
%              dimensionality over the batch.
%            h A NX1 cell array of function handles for the measurement
%              function that transform the state into the measurement
%              domain at each step. If the same measurement function is
%              used for all steps in the batch, then h can just be the
%              single function handle used.
%            F An xDim X xDim X (N-1) hypermatrix of  matrices. The state
%              at discrete-time k+1 is modeled as F(:,:,k) times the state
%              at time k plus zero-mean Gaussian process noise with
%              covariance matrix Q(:,:,k). Alternatively, if all of the
%              state transition matrices are the same, one can just pass a
%              single xDim X xDim matrix. Note that all of the F matrices
%              must be invertible if kD is not equal to one.
%            R The zDim X zDim X N hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a
%              single zDim X zDim matrix.
%           kD The discrete time-step at which the smoothed state estimate
%              is desired, where z(:,1) is at discrete time-step 1 (not 0).
%            Q The xDim X xDim X (N-1) hypermatrix of process noise
%              covariance matrices. Alternatively, if all of the process
%              noise covariance matrices are the same, one can just pass a
%              single xDim X xDim matrix. Note that all of the Q matrices
%              must be invertible except that zero matrices can be passed
%              (no process noise). If Q is omitted or an empty matrix is
%              passed, then the estimation is performed assuming there is
%              no process noise.
%       HJacob If a covariance matrix is desired on the output, then
%              HJacob can be provided. HJacob is a zDimX1 cell array of 
%              function handles for the measurement Jacobian matrix that
%              each takes the target state as a parameter. If a single
%              measurement Jacobian matrix is used for all steps of the
%              batch, then HJacob can just be the single function handle
%              used. If an empty matrix is passed or the parameter is
%              omitted, then HJacob will be found using numerical
%              differentiation via the numDiff function with default
%              parameters.
%   optimParam An optional structure whose members are arguments to the
%              function LSEstLMarquardt if one does not wish to use the
%              default values for that function. For example, to set
%              the maximum number of iterations to 1000, then use
%              optimParam.maxIter=1000; Possible members of the structure
%              are TolG,TolX,delta,deltaAbs,maxIter, and maxTries, which
%              are all described in the comments to the function
%              LSEstLMarquardt.
%
%OUTPUTS: xEst The refined batch state estimate at step kD. If the
%              optimization failed, then an empty matrix is returned.
%         PEst A covariance matrix estimate that goes with the state
%              estimate. The covariance matrix estimate is obtained using
%              the function batchLSLinMeasLinDyn. If the optimization
%              failed, then an empty matrix is returned.
%    xBatchEst The batch state estimates at all time-steps.
%     exitCode The exit status of the Levenburg-Marquart algorithm. The
%              values of this output are described in the comments to the
%              function LSEstLMarquardt.
%
%The algorithm is an implementation of the nonlinear ML technique discussed
%in  [1]. The aforementioned paper does not explicitly say how to perform
%the necessary least-squares optimization. Here, a form of the
%Levenburg-Marquart algorithm is used that does not require the evaluation
%of any gradients. If a covariance matrix estimate is desired, after a
%batch estimate has been obtained, the measurement function is linearized
%about the estimated state and the linear covariance estimation algorithm
%of batchLSLinMeasLinDyn is used.
%
%If process noise is included, then the estimation problem is over all of
%the states in the batch, though only the state at discrete time kD is
%returned. If process noise is omitted, then the estimation problem is only
%over the state at time kD.
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

    numSteps=size(z,2);
    xDim=size(xInit,1);
    zDim=size(z,1);
    
    if(nargin<6||isempty(kD))
        kD=1;
    end

    if(nargin<7)
        Q=[];
    end
    
    if(nargin<8)
        HJacob=[];
    end
    
    if(nargin<9)
        optimParam=[];
    end
    
    if(isa(h,'function_handle'))
       h=repmat({h},[numSteps,1]);
    end
    
    if(size(F,3)==1)
    	F=repmat(F,[1,1,numSteps-1]);
    end

    if(size(R,3)==1)
        SR=repmat(chol(R,'lower'),[1,1,numSteps]);
    else
        SR=zeros(size(R));
        for cStep=1:numSteps
            SR(:,:,cStep)=chol(R(:,:,cStep),'lower');
        end
    end
    
    %Get the values to pass to the LSEstLMarquardt function.
    if(~isempty(optimParam))
        if(isfield(optimParam,'TolG'))
            TolG=optimParam.TolG;
        else
            TolG=[];
        end

        if(isfield(optimParam,'TolX'))
            TolX=optimParam.TolX;
        else
            TolX=[];
        end

        if(isfield(optimParam,'delta'))
            delta=optimParam.delta;
        else
            delta=[];
        end

        if(isfield(optimParam,'deltaAbs'))
            deltaAbs=optimParam.deltaAbs;
        else
            deltaAbs=[];
        end

        if(isfield(optimParam,'maxIter'))
            maxIter=optimParam.maxIter;
        else
            maxIter=[];
        end

        if(isfield(optimParam,'maxTries'))
            maxTries=optimParam.maxTries;
        else
            maxTries=[];
        end
    end
    
    
    %If process noise is used, then estimation is performed over
    %the entire batch.
    if(~isempty(Q))
        if(size(Q,3)==1)
            %If Q is not all zeros
            if(~all(all(~Q)))
                SQ=repmat(chol(Q,'lower'),[1,1,numSteps-1]);
            else
                SQ=repmat(Q,[1,1,numSteps-1]);
            end
        else
            SQ=zeros(size(R));
            for cStep=1:numSteps
                if(~all(all(~Q(:,:,cStep))))
                    SQ(:,:,cStep)=chol(Q(:,:,cStep),'lower');
                else
                    SQ(:,:,cStep)=Q(:,:,cStep);
                end
            end
        end
        
        %If initial estimates over the entire batch were provided.
        if(size(xInit,2)>1)
            xBatchInit=xInit;
        else
            %The inital estimate over the batch is just xEst predicted forward and
            %backward in time.
            xBatchInit=zeros(xDim,numSteps);
            %The initial estimate.
            xBatchInit(:,kD)=xInit;
            %Predict forward in time
            for cStep=(kD+1):numSteps
                xBatchInit(:,cStep)=F(:,:,cStep-1)*xBatchInit(:,cStep-1);
            end

            %Predict backwards in time
            for cStep=(kD-1):-1:1
                xBatchInit(:,cStep)=F(:,:,cStep)\xBatchInit(:,cStep+1);
            end
        end
        
        if(isempty(optimParam))
            [xBatchEst,exitCode]=LSEstLMarquardt(@LMSCost,xBatchInit(:),false);
        else
            [xBatchEst,exitCode]=LSEstLMarquardt(@LMSCost,xBatchInit(:),false,TolG,TolX,delta,deltaAbs,maxIter,maxTries);
        end
        
        %If the optimization failed.
        if(isempty(xBatchEst))
           PEst=[];
           return;
        end
        
        xBatchEst=reshape(xBatchEst,xDim,numSteps);
        xEst=xBatchEst(:,kD);
    else
        %Otherwise, estimation is only over the single initial state.
        if(isempty(optimParam))
            [xEst,exitCode]=LSEstLMarquardt(@LMSCostNoProcNoise,xInit,false);
        else
            [xEst,exitCode]=LSEstLMarquardt(@LMSCostNoProcNoise,xInit,false,TolG,TolX,delta,deltaAbs,maxIter,maxTries);
        end
        
        %If the optimization failed.
        if(isempty(xEst))
           PEst=[];
           return;
        end
        
        %If estimates over the entire batch are needed for covariance
        %estimation
        if(nargout>2)
            xBatchEst=zeros(xDim,numSteps);
            %The refined initial estimate.
            xBatchEst(:,kD)=xEst;
            %Predict forward in time
            for cStep=(kD+1):numSteps
                xBatchEst(:,cStep)=F(:,:,cStep-1)*xBatchEst(:,cStep-1);
            end

            %Predict backwards in time
            for cStep=(kD-1):-1:1
                xBatchEst(:,cStep)=F(:,:,cStep)\xBatchEst(:,cStep+1);
            end
        end
    end

    if(nargout>1)
    %If a covariance matrix estimate is desired, use the linearized
    %estimate.
        if(nargin<8||isempty(HJacob))
            HJacob=cell(numSteps,1);
            for cStep=1:numSteps
                HJacob{cStep}=@(x)numDiff(x,h{cStep},zDim);
            end
        end

        if(isa(HJacob,'function_handle'))
            HJacob=repmat({HJacob},[numSteps,1]);
        end

        H=zeros(zDim,xDim,numSteps);
        for cStep=1:numSteps
            H(:,:,cStep)=HJacob{cStep}(xBatchEst(:,cStep));
        end

        [~,PEst]=batchLSLinMeasLinDyn([],H,F,R,kD,Q,numSteps);
    end
    
    function J=LMSCost(xBatch)
        %The cost function.
        
        xBatch=reshape(xBatch,xDim,numSteps);
        
        J=zeros(zDim*numSteps+xDim*(numSteps-1),1);
        
        diff=z(:,1)-h{1}(xBatch(:,1));
        J(1:zDim,1)=SR(:,:,1)\diff;
        
        minIdx=zDim+1;
        for curStep=2:numSteps
            diff=h{curStep}(xBatch(:,curStep))-z(:,curStep);
            
            maxIdx=minIdx+zDim-1;
            J(minIdx:maxIdx,1)=SR(:,:,curStep)\diff;
            minIdx=maxIdx+1;
            
            diff=xBatch(:,curStep)-F(:,:,curStep-1)*xBatch(:,curStep-1);
            maxIdx=minIdx+xDim-1;
            J(minIdx:maxIdx,1)=SQ(:,:,curStep-1)\diff;
            minIdx=maxIdx+1;
        end
    end

    function J=LMSCostNoProcNoise(xEst)
        J=zeros(zDim*numSteps,1);
        
        diff=h{kD}(xEst)-z(:,kD);
        J(1:zDim,1)=SR(:,:,kD)\diff;
        
        %Go forward in time from xEst
        xCur=xEst;
        minIdx=zDim+1;
        for curStep=(kD+1):numSteps
            xCur=F(:,:,curStep-1)*xCur;
            
            maxIdx=minIdx+zDim-1;
            diff=h{curStep}(xCur)-z(:,curStep);
            J(minIdx:maxIdx,1)=SR(:,:,curStep)\diff;
            minIdx=maxIdx+1;
        end
        
        %Go backward in time from xEst
        xCur=xEst;
        for curStep=(kD-1):-1:1
            xCur=F(:,:,curStep)\xCur;
            
            maxIdx=minIdx+zDim-1;
            diff=h{curStep}(xCur)-z(:,curStep);
            J(minIdx:maxIdx,1)=SR(:,:,curStep)\diff;
            minIdx=maxIdx+1;
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
