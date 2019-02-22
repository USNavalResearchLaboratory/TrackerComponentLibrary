function [xPred, PPred]=discEKFPred(xPrev,PPrev,f,FJacob,Q,FHessian)
%DISCEKFPRED Perform the discrete-time prediction step that comes with 
%            the first- or second-order order extended Kalman filter (EKF).
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        PPrev The xDim X xDim state covariance matrix at the previous
%              time-step.
%            f A function handle for the state transition function that
%              takes the state as its parameter.
%       FJacob A function handle for calculating the xDim X xDim Jacobian
%              of f, or the xDim X xDim Jacobian matrix itself. If an empty
%              matrix is passed, then FJacob will be found using numerical
%              differentiation  via the numDiff function with default
%              parameters.
%            Q The xDimX xDim process noise covariance matrix.
%     FHessian This parameter is only provided if a second-order EKF is
%              desired. This is either a function handle for the state
%              transition Hessian hypermatrix, or it is the state
%              transition Hessian hypermatrix itself. The matrix is
%              xDim X xDim X xDim. The matrix FH=FHessian(x) is such that
%              FH(i,j,k) is the second derivative of the kth element of the
%              vector returned by f with respect to the ith and jth
%              components of x. The Hessian matrix is symmetric. If this
%              parameter is omitted, a first-order EKF is used.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate.
%         PPred The xDim X xDim predicted state covariance matrix.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of
%[1]. The second-order EKF is described in Chapter 10.3.2 of the same text.
%
%The partial derivatives in the Jacobain matrix returned by the function
%FJacob are ordered
%[dF/dx(1), dF/dx(2),...,dF/dx(xDim)]
%That is, column i consists of partial derivatives with respect to element
%i of the x vector.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
if(isempty(FJacob))
    FJacob=@(x)numDiff(x,f,xDim);
end

if(nargin<6)
   FHessian=[]; 
end

if(isa(FJacob,'function_handle'))
    F=FJacob(xPrev);
else
    F=FJacob;%If the Jacobian matrix was directly given.
end

if(~isempty(FHessian))
    if(isa(FHessian,'function_handle'))
        FHess=FHessian(xPrev);
    else
        FHess=FHessian;%If the Hessian matrix was directly given.
    end

    %If terms using the Hessian (for the second-order filter) should be
    %computed.
    PxxHessTerm=zeros(xDim,xDim);
    xPredHessTerm=zeros(xDim,1);
    for n=1:xDim
        en=zeros(xDim,1);
        en(n)=1;

        HPProdn=FHess(:,:,n)*PPrev;
        xPredHessTerm=xPredHessTerm+en*trace(HPProdn);

        for m=1:xDim
            em=zeros(xDim,1);
            em(m)=1;
            PxxHessTerm=PxxHessTerm+en*em'*trace(HPProdn*FHess(:,:,m)'*PPrev);
        end
    end
    PxxHessTerm=PxxHessTerm/2;
    xPredHessTerm=xPredHessTerm/2;
else
    PxxHessTerm=0;
    xPredHessTerm=0;
end
xPred=f(xPrev)+xPredHessTerm;
PPred=F*PPrev*F'+Q+PxxHessTerm;
%Handle possible loss of symmetry due to order of operations and finite
%precision limitations.
PPred=(PPred+PPred')/2;
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
