function [xEst,PEst]=batchLSNonlinMeasNonlinDyn(xInit,z,h,R,HJacob,numIter)
%%BATCHLSNONLINMEASNONLINDYN Perform batch least squares state estimation
%                      refinement under a nonlinear measurement model and a
%                      linear dynamic model without process noise. The
%                      iterative nonlinear least squares algorithm is used.
%
%INPUTS: xInit The initial estimate of the target state at the time of the
%              kDth measurement that is to be refined through iteration.
%            z The zDim X N matrix of measurements for the whole batch. It
%              is assumed that the measurements have the same
%              dimensionality over the batch.
%            h A NX1 cell array of function handles for the measurement
%              function that transform the state into the measurement
%              domain at each step. If the same measurement function is
%              used for all steps in the batch, then h can just be the
%              single function handle used.
%            R The zDim X zDim X N hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a
%              single zDim X zDim matrix.
%       HJacob A NX1 cell array of function handles for the measurement
%              Jacobian matrix that each takes the target state as a
%              parameter. If a single measurement Jacobian matrix is used
%              for all steps of the batch, then HJacob can just be the
%              single function handle used. If an empty matrix is passed or
%              HJacob is not provided, then HJacob will be found using
%              numerical differentiation via the numDiff function with
%              default parameters.
%      numIter The numer of iterations to perform in the nonlinear least
%              squares algorithm. If this parameter is omitted, the
%              default is 10.
%
%OUTPUTS: xEst The refined batch state estimate at step kD.
%         PEst A covariance matrix estimate that goes with the state
%              estimate.
%
%The algorithm is an implementation of the nonlinear iterative least
%squares algorithm of Chapter 3.4 of [1]. The covariance provided is the
%covariance associated with the algorithm. The measurement matrices for the
%propagated state as well as the transition matrices for the propagated
%state are taken as derived in Section 3.3.2 of [2].
%
%The algorithm is likely to diverge if xEst is particularly bad. A more
%robust albeit slower estimation algorithm would use a Levenburg-Marquart
%algorithm to assure that each step improves the cost function.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] A. B. Poore, B. J. Slocumb, B. J. Suchomel, F. H. Obermeyer, S. M.
%    Herman, and S. M. Gadaleta, "Batch maximum likelihood (ML) and maximum
%    a posteriori (MAP) estimation with process noise for tracking
%    applications," in Proceedings of SPIE: Signal and Data Processing of
%    Small Targets, vol. 5204, San Diego, CA, 3 Aug. 2003, pp. 188-199.
%
%January 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(numIter))
    numIter=10;
end
if(nargin<5)
    HJacob=[];
end

xEst=xInit;

zDim=size(z,1);
xDim=size(xEst,1);
numSteps=size(z,2);

if(isa(h,'function_handle'))
    h=repmat({h},[numSteps,1]);
end

if(isa(HJacob,'function_handle'))
    HJacob=repmat({HJacob},[numSteps,1]);
end

if(isempty(HJacob))
    HJacob=cell(numSteps,1);
    for curStep=1:numSteps
        HJacob{curStep}=@(x)numDiff(x,h{curStep},zDim);
    end
end

if(size(R,3)==1)
    R=repmat(R,[1,1,numSteps]);
end

J=zeros(zDim*numSteps,xDim);
zPredStacked=zeros(zDim*numSteps,1);

RStackedInv=inv(blkDiagRep(R));
for curIter=1:numIter
    minIdx=1;
    for curStep=1:numSteps
        maxIdx=minIdx+zDim-1;
        
        J(minIdx:maxIdx,:)=HJacob{curStep}(xEst);
        
        zPredStacked(minIdx:maxIdx)=h{curStep}(xEst);
        minIdx=maxIdx+1;
    end
    
    xEst=xEst+lsqminnorm(J'*RStackedInv*J,J'*RStackedInv*(z(:)-zPredStacked));
end

%If a covariance matrix estimate is desired.
if(nargout==2)
    minIdx=1;
    for curStep=1:numSteps
        maxIdx=minIdx+zDim-1;
        
        J(minIdx:maxIdx,:)=HJacob{curStep}(xEst);
        
        minIdx=maxIdx+1;
    end
    
    PEst=pinv(J'*RStackedInv*J);
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
