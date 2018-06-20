function [xEst,PEst,exitCode]=batchLSNonlinMeasNonlinDynLM(xInit,z,h,R,kD,HJacob,optimParam)
%%BATCHLSNONLINMEASLINDYNLM  Given a single measurement or a batch of
%           measurements, perform batch least squares state estimation
%           refinement under a nonlinear measurement model and a linear
%           dynamic model. Nonlinear least squares optimization is 
%           performed using the Levenburg-Marquart algorithm with a 
%           numerically estimated Jacobian.
%
%INPUTS: xInit The xDimX1 initial estimate of the target state at the time
%              of the kDth measurement.
%            z The zDim X N matrix of measurements for the whole batch. It
%              is assumed that the measurements have the same
%              dimensionality over the batch.
%            h A NX1 cell array of function handles for the measurement
%              function that transform the state into the measurement
%              domain at each step. If the same measurement function is
%              used for all steps in the batch, then h can just be the
%              single function handle used. It is assumed that any 
%              required state propagation forward and backward from time-
%              step kd is included.
%            R The zDim X zDim X N hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a
%              single zDim X zDim matrix.
%           kD The discrete time-step at which the smoothed state estimate
%              is desired, where z(:,1) is at discrete time-step 1 (not
%              0). If kD is omitted or an empty matrix is passed, a value
%              of 1 is assumed.
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
%              estimate. If the optimization failed, then an empty matrix
%              is returned.
%     exitCode The exit status of the Levenburg-Marquart algorithm. The
%              values of this output are described in the comments to the
%              function LSEstLMarquardt.
%
%The algorithm is an implementation of the nonlinear ML technique discussed
%in [1]. The aforementioned paper does not explicitly say how to perform
%the necessary least-squares optimization. Here, a form of the
%Levenburg-Marquart algorithm from immoptitoolbox is used that does not
%require the evaluation of any gradients. If a covariance matrix estimate
%is desired, after a batch estimate has been obtained, the measurement
%function is linearized about the estimated state and the linear covariance
%estimation algorithm of batchLSLinMeasLinDyn is used.
%
%REFERENCES:
%[1] A. B. Poore, B. J. Slocumb, B. J. Suchomel, F. H. Obermeyer, S. M.
%    Herman, and S. M. Gadaleta, "Batch maximum likelihood (ML) and maximum
%    a posteriori (MAP) estimation with process noise for tracking
%    applications," in Proceedings of SPIE: Signal and Data Processing of
%    Small Targets, vol. 5204, San Diego, CA, 3 Aug. 2003, pp. 188-199.
%
%Febraury 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(kD))
    kD=1;
end

if(nargin<6)
    optimParam=[];
end

numSteps=size(z,2);
zDim=size(z,1);

if(isa(h,'function_handle'))
    h=repmat({h},[numSteps,1]);
end

if(size(R,3)==1)
    SR=repmat(chol(R,'lower'),[1,1,numSteps]);
else
    SR=zeros(size(R));
    for cStep=1:numSteps
        SR(:,:,cStep)=chol(R(:,:,cStep),'lower');
    end
end

if(isempty(optimParam))
    [xEst,exitCode]=LSEstLMarquardt(@LMSCostNoProcNoise,xInit);
else
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
    
    [xEst,exitCode]=LSEstLMarquardt(@LMSCostNoProcNoise,xInit,false,TolG,TolX,delta,deltaAbs,maxIter,maxTries);
end
    
%If the optimization failed.
if(isempty(xEst))
   PEst=[];
   return;
end

%If the covariance is being calculated
if(nargout>1)
    if(nargin<6||isempty(HJacob))
        HJacob=cell(numSteps,1);
        for cStep=1:numSteps
            HJacob{cStep}=@(x)numDiff(x,h{cStep},zDim);
        end
    end
    if(isa(HJacob,'function_handle'))
        HJacob=repmat({HJacob},[numSteps,1]);
    end
    
    RStackedInv=inv(blkDiagRep(SR));
    minInd=1;
    for cStep=1:numSteps
        maxInd=minInd+zDim-1;
        
        J(minInd:maxInd,:)=HJacob{cStep}(xEst);
        
        minInd=maxInd+1;
    end
    
    PEst=pinv(J'*RStackedInv*J);
end

    function J=LMSCostNoProcNoise(xEst)
        J=zeros(zDim*numSteps,1);
        
        diff=h{kD}(xEst)-z(:,kD);
        J(1:zDim,1)=SR(:,:,kD)\diff;
        
        %Go forward in time from xEst
        minIdx=zDim+1;
        for curStep=(kD+1):numSteps
            maxIdx=minIdx+zDim-1;
            diff=h{curStep}(xEst)-z(:,curStep);
            J(minIdx:maxIdx,1)=SR(:,:,curStep)\diff;
            minIdx=maxIdx+1;
        end
        
        %Go backward in time from xEst
        for curStep=(kD-1):-1:1
            maxIdx=minIdx+zDim-1;
            diff=h{curStep}(xEst)-z(:,curStep);
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
