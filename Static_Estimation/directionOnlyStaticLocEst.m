function [t,exitCode]=directionOnlyStaticLocEst(u,lRx,algorithm,params1,params2,params3)
%%DIRECTIONONLYSTATICLOCEST Given simultaneous direction estimates towards
%                      a target from at least two sensors, estimate the
%                      Cartesian location of the target in 2D or 3D.
%
%INPUTS: u A 2XnumMeas or 3XnumMeas set of 2D or 3D unit direction vectors 
%          in a global coordinate system pointing from each sensor to the
%          target in 2D or 3D. numMeas>=2.
%      lRx The 2XnumMeas or 3XnumMeas set of Cartesian locations of the
%          sensors producing the direction vector measurements.
% algorithm A parameter that specifies how the location of the target is to
%          be estimated. The algorithm chosen specifies what the
%          parameters params1, params2, and params3 are. Possible values
%          are:
%          0 (The default if omitted or an empty ,matrix is passed) Use
%            the suboptimal least-squares triangulation algorithm. This
%            may be followed by iterations of an explicit solution to the
%            maximum likelihood conditioned on the distances from each
%            receiver to the target being known.
%          1 Use the suboptimal least-squares triangulation algorithm to
%            get an initial estimate, which is provided to a quasi-Newton
%            method to locally maximize the likelihood.
%          2 Use an explicit solution given the distances from each
%            receiver to the target.
%          3 Use a quasi-Newton method to maximize the likelihood given an
%            initial estimate.
%  params1 The first set of parameters for the selected algorithm. The
%          meaning depends on the value of algorithm. params1 has the
%          following meaning for the following values of algorithm:
%           algorithm=0,1: params1 is the parameters for the suboptimal
%             least-squares triangulation algorithm. params1 is a
%             structure that can contain elements with the following names
%             and meanings:
%             'W' If weight matrices are to be used in the suboptimal
%                 least squares algorithm (not recommended), then W is a
%                 2X2XnumMeas (for 2D) or 3X3XnumMeas (for 3D) set of the
%                 weights for use in the cost function. If this parameter
%                 is omitted or an empty matrix is passed, the unweighted
%                 algorithm is used (W is all identity matrices).
%             useConstAlg A boolean parameter indicating whether a
%                 constrained version of the suboptimal least squares
%                 algorithm should be used. Generally, the simpler,
%                 unconstrained algorithm suffices. The constraints make
%                 sure that the estimate is chosen such that all range
%                 estimates used in the algorithm from the sensors to the
%                 target are positive. The default if omitted or an empty
%                 matrix is passed is false.
%           algorithm=2: In this instance, params1 is the parameters for
%                 the explicit solution given the ranges from the sensors
%                 to the target. The element 'r' of the structure is not
%                 optional. Possible elements are
%                 'r' A numMeasX1 or 1XnumMeas vector of the ranges from
%                     the target to each sensor (for each measurement).
%                 'RInv' A 2X2XnumMeas (for 2D) or 3X3XnumMeas (for 3D)
%                     set of inverse covariance matrices of the
%                     measurements. The matrices can be singular. If this
%                     parameter is omitted or an empty matrix is passed,
%                     then identity matrices shall be used.
%                 'numIter' The number of iterations of the algorithm to
%                     perform. Each iteration after the first is the same
%                     algorithm, except it is initialized using the range
%                     estimates computed from the previous estimate. The
%                     default if this parameter is omitted or an empty
%                     matrix is passed is 1. numIter>=0.
%           algorithm=3: In this instance, params1 is a structure of the
%                 parameters for a quasi-Newton method to maximize the
%                 likelihood. The initial estimate 'tInit' is required.
%                 Possible values are:
%                 'tInit' The 2X1 or 3X1 initial estimate of the target
%                     location.
%                 'RInv' A 2X2XnumMeas (for 2D) or 3X3XnumMeas (for 3D)
%                     set of inverse covariance matrices of the
%                     measurements. The matrices can be singular. If this
%                     parameter is omitted or an empty matrix is passed,
%                     then identity matrices shall be used. 
%                 'epsilon','deltaTestDist','delta','lineSearchParams',
%                  'scaleD','maxIter' -All of these parameters correspond
%                     to the same named inputs of the function
%                     quasiNetwonBFGS. See the comments to the function
%                     quasiNetwonBFGS for more details.
%  params2 The second set of parameters for the selected algorithm. The
%          meaning depends on the value of algorithm. params2 has the
%          following meaning for the following values of algorithm:
%          algorithm=0 This is the same as params1 for algorithm=2, except
%                  no initial estimate 'r' will be used as the one
%                  provided by the suboptimal least squares algorithm will
%                  be used.
%          algorithm=1 This is the same as params1 for algorithm=2, except
%                  'tInit' is not used as the initial target location
%                  estimate will be obtained using the suboptimal least
%                  squares algorithm.
%          algorithm=2,3 params2 is not used, does not need to be
%                  provided, and will be ignored if provided.
%  params3 These parameters are only used if algorithms 0 or 1 are chosen
%          AND params1.useConstAlg is true. params3 is a structure holding
%          the parameters for the function convexQuadProg. Possible
%          entries are 'epsVal', and 'maxIter' and are as described in the
%          comments to convexQuadProg.
%
%OUTPUTS: t The 2X1 or 3X1 estimate of the target location.
%  exitCode A 2X1 vector indicating how the algorithm terminated. When
%           exitCode(1)=0 means that the second half indicates that the
%           value of exitCode(2) represents a value returned by the
%           convexQuadProg function (or 0 is no error) and exitCode(1)=1
%           indicates that exitCode(2) represents a value returned by the
%           quasiNetwonBFGS function.
%
%The algorithms are based on [1].
%
%The quasi-Newton method estimator minimizes the cost function
%\sum_{i=1}^numMeas((t-lRx(:,i))/norm((t-lRx(:,i)))-u(:,i))'*RInv(:,:,i)*((t-lRx(:,i))/norm((t-lRx(:,i)))-u(:,i))
%over t using the quasiNetwonBFGS algorithm.
%
%The explicit solution conditioned on knowing the distances from the target
%to the sensors that minimizes the cost function
%\sum_{i=1}^numMeas((t-lRx(:,i))/r(i)-u(:,i))'*RInv(:,:,i)*((t-lRx(:,i))/r(i)-u(:,i))
%over t where r(i) is the known target-to-sensor distance.
%
%The suboptimal least squares solution minimizes the cost function
%\sum_{i=1}^numMeas(t-lRx(:,i)-r(i)*u(:,i))'*RInv(:,:,i)*(t-lRx(:,i)-r(i)*u(:,i))
%where the optimization is over t AND over all r(i). If the constrained
%algorithm is used, then the condition that all r(i)>=0 is enforced. The
%constrained algorithm uses the convexQuadProg function.
%
%EXAMPLE 1:
%As an example, consider 4 sensors  and one target. The sensor take
%measurements in a local u-v coordinate system The results must be
%augmented to become full unit vectors and rotated into the global
%coordinate system.
%Sensors on Hawaiian islands
% lRxLatLonAlt=[[20.72529087;-156.1164093;1618],...
%               [20.74070316;-155.98731995;36],...
%               [20.25189031;-155.80604553;53],...
%               [20.1152606;-155.54855347;187]];
% lRxLatLonAlt(1:2,:)=lRxLatLonAlt(1:2,:)*(pi/180);%Convert to radians.
% lRx=ellips2Cart(lRxLatLonAlt);%Convert to Cartesian
% targetLatLonAlt=[20.62250226;-155.42495728;8000];%Middle scenario
% targetLatLonAlt(1:2)=targetLatLonAlt(1:2)*(pi/180);%Convert to radians.
% targetLoc=ellips2Cart(targetLatLonAlt);%Convert to Cartesian
% %Create the unit vector from each receiver to the target.
% u=bsxfun(@minus,targetLoc,lRx);
% u=bsxfun(@rdivide,u,sqrt(sum(u.*u,1)));
% 
% %Two radars pointed East, Two radars pointed North.
% az12=90*(pi/180);%East
% el=15*(pi/180);%Radar elevation, radians.
% %M will be rotation matrices to go from global to local coordinates.
% M=zeros(3,3,4);
% M(:,:,1)=findRFTransParam(lRxLatLonAlt(:,1),az12,el);
% M(:,:,2)=findRFTransParam(lRxLatLonAlt(:,2),az12,el);
% %The final two radars are pointing North
% az34=0*(pi/180);%North
% M(:,:,3)=findRFTransParam(lRxLatLonAlt(:,3),az34,el);
% M(:,:,4)=findRFTransParam(lRxLatLonAlt(:,4),az34,el);
% 
% numSensors=4;
% lRx=lRx(:,1:numSensors);
% %The noise-free target directions in the local coordinate systems of the
% %receivers are
% uRx=zeros(3,numSensors);
% for i=1:numSensors
%     %M to go from global to local
%     uRx(:,i)=M(:,:,i)*u(:,i);
% end
% 
% SR=diag([1e-3;1e-3]);%Square-root Measurement covariance matrix in local
%                      %u-v coordinates. The error is about 1mrad in azimuth
%                      %and elevation when near the boresight.
% R=SR*SR';%The covariance matrix.
% RInv=inv(R);  
% %Get the inverse covariance matrix in global coordinates.
% RInvGlobal=zeros(3,3,4);%Allocate space
% for i=1:numSensors
%     %The third coordinate provides no additional information locally, so
%     %the inverse covariance values for it are set to zero (a singular
%     %inverse covariance matrix).
%     RInvGlobal(1:2,1:2,i)=RInv;
%     RInvGlobal(:,:,i)=M(:,:,i)*RInvGlobal(:,:,i)*M(:,:,i)';
% end
% 
% %Measurements
% uMeas=zeros(3,numSensors);
% uMeas(1:2,:)=uRx(1:2,:)+SR*randn(2,numSensors);%Local measurement.
% uMeas(3,:)=sqrt(1-sum(uMeas(1:2,:).^2,1));%Make a unit vector in 3D
% for i=1:numSensors
%     uMeas(:,i)=M(:,:,i)'*uMeas(:,i);%Rotate into global coordinates.
% end
% [t,exitCode]=directionOnlyStaticLocEst(uMeas,lRx,0);
% estError=norm(t-targetLoc)%The estimation error.
%
%EXAMPLE 2:
%When this function is called with params2=[];, params2.numIter=0; the
%output is the same as if one were to call closestPointBetween2Lines to
%just find the closest point between lines going out from each receiver. We
%show this by showing that the difference between the results is zero.
% l1=[0;0;0];
% l2=[100;0;0];
% tarLoc=[200;100;50];
% u1=(tarLoc-l1)./norm(tarLoc-l1);
% u2=(tarLoc-l2)./norm(tarLoc-l2);
% %Add some noise so there isn't a perfect intersection.
% u1=u1+randn(3,1);
% u1=u1/norm(u1);
% u2=u2+randn(3,1);
% u2=u2/norm(u2);
% 
% params2=[];
% params2.numIter=0;
% tLocEst=directionOnlyStaticLocEst([u1,u2],[l1,l2],0,[],params2);
% pointBetweenLines=closestPointBetween2Lines(u1,l1,u2,l2);
% AbsDiff=tLocEst-pointBetweenLines
%
%REFERENCES:
%[1] D. F. Crouse, "Bearings-Only Localization Using Direction Cosines," in
%    Proceedings of the 19th International Conference on Information
%    Fusion, Heidelberg, Germany, July 2016.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

numDim=size(u,1);
numMeas=size(u,2);

exitCode=[0;0];

switch(algorithm)
    case 0%Suboptimal least squares followed by iterations of the explicit
          %algorithm.

        %Default parameters if params1 does not supply them.
        W=[];
        useConstAlg=false;
          
        if(nargin>3&&~isempty(params1))
            %Extract any parameters passed.
            if(isfield(params1,'W'))
               W=params1.W; 
            end
            
            if(isfield(params1,'useConstAlg'))
               useConstAlg=params1.useConstAlg; 
            end
        end
        
        %Default parameters if params2 does not supply them.
        RInv=repmat(ones(numDim,numDim),1,1,numMeas);
        numIter=1;
        
        if(nargin>4&&~isempty(params2))
            if(isfield(params2,'RInv'))
               RInv=params2.RInv; 
            end
            
            if(isfield(params2,'numIter'))
               numIter=params2.numIter; 
            end
        end
        
        if(nargin<6)
            params3=[];
        end
        
        %Call the suboptimal least-squares algorithm.
        [t,retVal]=suboptimalLSTriangulation(u,lRx,W,useConstAlg,params3);
        if(retVal~=0)%If an error occurred.
            exitCode(2)=retVal;
            return;
        end

        %Refine the estimate, if possible, using the explicit algorithm.
        for curIter=1:numIter
            r=sqrt(sum(bsxfun(@minus,t,lRx).^2,1));
            t=triangulateKnownR(r,u,lRx,RInv);
        end
    case 1%Suboptimal least squares followed by iterations of the quasi-
          %Newton method.
          
        %Default parameters if params1 does not supply them.
        W=[];
        useConstAlg=false;
          
        if(nargin>3&&~isempty(params1))
            %Extract any parameters passed.
            if(isfield(params1,'W'))
               W=params1.W; 
            end
            
            if(isfield(params1,'useConstAlg'))
               useConstAlg=params1.useConstAlg; 
            end
        end
        
        %Default parameters if params2 does not supply them.
        RInv=repmat(ones(numDim,numDim),1,1,numMeas);
        epsilon=[];
        deltaTestDist=[];
        delta=[];
        lineSearchParams=[];
        scaleD=[];
        maxIter=[];
        
        if(nargin>4&&~isempty(params2))
            if(isfield(params2,'RInv'))
               RInv=params2.RInv; 
            end
            
            if(isfield(params2,'epsilon'))
               epsilon=params2.epsilon; 
            end
            
            if(isfield(params2,'deltaTestDist'))
               deltaTestDist=params2.deltaTestDist; 
            end
            
            if(isfield(params2,'delta'))
               delta=params2.delta; 
            end
            
            if(isfield(params2,'lineSearchParams'))
               lineSearchParams=params2.lineSearchParams; 
            end
            
            if(isfield(params2,'scaleD'))
               scaleD=params2.scaleD; 
            end
            
            if(isfield(params2,'maxIter'))
               maxIter=params2.maxIter; 
            end
        end
        
        if(nargin<6)
            params3=[];
        end
        
        %Call the suboptimal least-squares algorithm.
        [t,retVal]=suboptimalLSTriangulation(u,lRx,W,useConstAlg,params3);
        
        if(retVal~=0)%If an error occurred.
            exitCode(2)=retVal;
            return;
        end
                
        %Call Newton's method using the least squares initial estimate.
        [t,retVal]=triangulateQuasiNewton(t,u,lRx,RInv,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter);
        if(retVal~=0)%If an error occurred.
            exitCode(1)=1;
            exitCode(2)=retVal;
            return;
        end
    case 2%The explicit algorithm requiring range information.
        %Default parameters if params1 does not supply them.
        RInv=repmat(ones(numDim,numDim),1,1,numMeas);
        numIter=1;
        
        r=params1.r;
        if(isfield(params1,'RInv'))
           RInv=params1.RInv; 
        end

        if(isfield(params1,'numIter'))
           numIter=params1.numIter; 
        end

        t=triangulateKnownR(r,u,lRx,RInv);
        %Refine the estimate, if possible, using the explicit
        for curIter=1:numIter
            t=triangulateKnownR(r,u,lRx,RInv);
            r=sqrt(sum(bsxfun(@minus,r,lRx).^2,1));
        end
    case 3%The quasi-Newton method.
        
        %Default parameters if params1 does not supply them.
        RInv=repmat(ones(numDim,numDim),1,1,numMeas);
        epsilon=[];
        deltaTestDist=[];
        delta=[];
        lineSearchParams=[];
        scaleD=[];
        maxIter=[];
        
        tInit=params1.tInit;
        if(isfield(params1,'RInv'))
           RInv=params1.RInv; 
        end

        if(isfield(params1,'epsilon'))
           epsilon=params1.epsilon; 
        end

        if(isfield(params1,'deltaTestDist'))
           deltaTestDist=params1.deltaTestDist; 
        end

        if(isfield(params1,'delta'))
           delta=params1.delta; 
        end

        if(isfield(params1,'lineSearchParams'))
           lineSearchParams=params1.lineSearchParams; 
        end

        if(isfield(params1,'scaleD'))
           scaleD=params1.scaleD; 
        end

        if(isfield(params1,'maxIter'))
           maxIter=params1.maxIter; 
        end
        
        %Call Newton's method using the provided initial estimate.
        [t,retVal]=triangulateQuasiNewton(tInit,u,lRx,RInv,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter);
        if(retVal~=0)%If an error occurred.
            exitCode(1)=1;
            exitCode(2)=retVal;
            return;
        end
    otherwise
        error('Unknown algorithm specified')
end


end

function [t,exitCode]=suboptimalLSTriangulation(u,lRx,W,useConstAlg,quadProgParams)
%%SUBOPTIMALLSTRIANGULATION This function implements the suboptimal least
%squares triangulation where the optimization uses the cost function
%\sum_{i=1}^numMeas(t-lRx(:,i)-r(i)*u(:,i))'*RInv(:,:,i)*(t-lRx(:,i)-r(i)*u(:,i))
%where the optimization is over t AND over all r(i).If the constrained
%algorithm is used, then the condition that all r(i)>=0 is enforced. The
%constrained algorithm uses the convexQuadProg function.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(useConstAlg))
    useConstAlg=false;
end

exitCode=0;

%The number of observations.
n=size(u,2);

%Use the unconstrained algorithm for W being all identity matrices if W is
%not provided or an empty matrix is passed for W.
if(nargin<3||isempty(W)&&useConstAlg==false)
    %The unconstrained solution.
    A=zeros(n,n);
    b=zeros(n,1);

    for i=1:n
        A(i,i)=1-1/n;
        for j=(i+1):n
            A(i,j)=-(1/n)*dot(u(:,i),u(:,j));
            A(j,i)=A(i,j);
        end
        b(i)=(1/n)*sum(sum(bsxfun(@times,lRx,u(:,i)),1))-dot(lRx(:,i),u(:,i));
    end

    r=A\b;
    
    t=(1/n)*sum(lRx+bsxfun(@times,r',u),2);
    return;
end

numDim=size(u,1);

%If W is not provided, then let it be all identity matrices
if(isempty(W))
   W=repmat(eye(numDim,numDim),1,1,n);
end

exitCode=0;

if(useConstAlg==false)
    %If the unconstrained algorithm should be used.
    A=zeros(n,n);
    b=zeros(n,1);
    
    WSumInv=inv(sum(W,3));
    WLSum=zeros(numDim,1);
    for i=1:n
        WLSum=WLSum+W(:,:,i)*lRx(:,i);
    end
    
    for i=1:n
        A(i,i)=u(:,i)'*W(:,:,i)*u(:,i)-u(:,i)'*W(:,:,i)*WSumInv*W(:,:,i)*u(:,i);
        for j=(i+1):n
            A(i,j)=-u(:,i)'*W(:,:,i)*WSumInv*W(:,:,j)*u(:,j);
            A(j,i)=A(i,j);
        end
        b(i)=u(:,i)'*W(:,:,i)*WSumInv*WLSum-u(:,i)'*W(:,:,i)*lRx(:,i);
    end
    r=A\b;
    
    t=zeros(numDim,1);
    for i=1:n
        t=t+W(:,:,i)*(lRx(:,i)+r(i)*u(:,i));
    end
    t=WSumInv*t;
    return;
else
    %Default parameters for the convexQuadProg function.
    epsVal=[];
    maxIter=[];
    
    if(~isempty(quadProgParams))
        if(isfield(quadProgParams,'epsVal'))
           epsVal=quadProgParams.epsVal; 
        end
        
        if(isfield(quadProgParams,'maxIter'))
           maxIter=quadProgParams.maxIter; 
        end
    end
    
    %Formulate the problem as a quadratic programming problem.
    Q=zeros(numDim+n,numDim+n);
    c=zeros(numDim+n,1);
    
    Q(1:numDim,1:numDim)=sum(W,3);
    WLSum=zeros(numDim,1);
    for i=1:n
        WLSum=WLSum+W(:,:,i)*lRx(:,i);
    end
    c(1:numDim,1)=-WLSum;
    
    for i=1:n
        Q(1:numDim,numDim+i)=-W(:,:,i)*u(:,i);
        Q(numDim+i,1:numDim)=Q(1:numDim,numDim+i)';
        Q(numDim+i,numDim+i)=u(:,i)'*W(:,:,i)*u(:,i);

        c(numDim+i)=lRx(:,i)'*W(:,:,i)*u(:,i);
    end
    
    C=[zeros(n,numDim),eye(n,n)]';
    b=zeros(n,1);
    [tTilde,~,exitCode]=convexQuadProg(Q,c,C,b,0,epsVal,maxIter);
    t=tTilde(1:numDim);
    %r=tTilde((numDim+1):end);
    return;
end

end

function t=triangulateKnownR(r,u,lRx,RInv)
%%TRIANGULATEKNOWNR This function implements the explicit solution
%conditioned on knowing the distances from the target to the sensors
%that minimizes the cost function
%\sum_{i=1}^numMeas((t-lRx(:,i))/r(i)-u(:,i))'*RInv(:,:,i)*((t-lRx(:,i))/r(i)-u(:,i))
%over t where r(i) is the known target-to-sensor distance.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(u,1);
numMeas=size(u,2);

t=zeros(numDim,1);

RInvSum=zeros(numDim,numDim);
for i=1:numMeas
    RInv(:,:,i)=eye(3);
    RInvSum=RInvSum+(1/r(i)^2)*RInv(:,:,i);
    
    t=t+(1/r(i))*RInv(:,:,i)*((1/r(i))*lRx(:,i)+u(:,i));
end
t=RInvSum\t;

end

function [t,exitCode]=triangulateQuasiNewton(tInit,u,lRx,RInv,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter)
%%TRIANGULATEQUASINEWTON This formulates the maximum likelihood cost
%function. The cost function
%\sum_{i=1}^numMeas((t-lRx(:,i))/norm((t-lRx(:,i)))-u(:,i))'*RInv(:,:,i)*((t-lRx(:,i))/norm((t-lRx(:,i)))-u(:,i))
%is to be minimized over t. An initial estimate tInit must be provided.
%The cost function used by this function is called costFunc and is in this
%same file.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

f=@(t)costFunc(t,u,lRx,RInv);

[t,~,exitCode]=quasiNetwonBFGS(f,tInit,[],epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter);
end

function [val,grad]=costFunc(t,u,lRx,RInv)
    numVals=size(u,2);

    val=0;
    for i=1:numVals
        tlDiff=t-lRx(:,i);
        ltDiffMag=norm(tlDiff);
        
        diff=tlDiff/ltDiffMag-u(:,i);
        
        val=val+diff'*RInv(:,:,i)*diff;
    end
    
    numDim=size(u,1);
    
    %For the gradient.
    grad=zeros(numDim,1);
    
    for i=1:numVals
        tlDiff=t-lRx(:,i);
        ltDiffMag=norm(tlDiff);
        
        if(numDim==2)
            A=[tlDiff(2)^2,            -tlDiff(1)*tlDiff(2);
               -tlDiff(1)*tlDiff(2),     tlDiff(1)^2];
        else
            A=[tlDiff(2)^2+tlDiff(3)^2, -tlDiff(1)*tlDiff(2),   -tlDiff(1)*tlDiff(3);
               -tlDiff(1)*tlDiff(2),    tlDiff(1)^2+tlDiff(3)^2,-tlDiff(2)*tlDiff(3);
               -tlDiff(1)*tlDiff(3),    -tlDiff(2)*tlDiff(3),   tlDiff(1)^2+tlDiff(2)^2];
        end
        
        grad=grad+(1/ltDiffMag^4)*A*RInv(:,:,i)*tlDiff-(1/ltDiffMag^3)*A*RInv(:,:,i)*u(:,i);
    end
    
    grad=2*grad;
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
