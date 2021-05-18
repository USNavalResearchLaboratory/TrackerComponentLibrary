function [xPred, SPred]=sqrtDiscEKFPred(xPrev,SPrev,f,FJacob,SQ)
%SQRTDISCEKFPRED Perform the discrete-time prediction step that comes with 
%                the square-root implementation of the first-order Extended
%                Kalman Filter (EKF).
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        SPrev The xDim X xDim lower-triangular square root of the  state
%              covariance matrix at the previous time-step.
%            f A function handle for the state transition function that
%              takes the state as its parameter.
%       FJacob A function handle for calculating the xDim X xDim state
%              transition matrix. If an empty matrix is passed, then FJacob
%              will be found using numerical differentiation  via the
%              numDiff function with default parameters.
%           SQ The xDimX xDim lower-triangular square root of the  process
%              noise covariance matrix.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate.
%         SPred The xDim X xDim lower-triangular square root of the
%               predicted state covariance estimate.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of
%[1].
%
%The partial derivatives in the Jacobian matrix returned by the function
%FJacob are ordered
%[dF/dx(1), dF/dx(2),...,dF/dx(xDim)]
%That is, column i consists of partial derivatives with respect to element
%i of the x vector.
%
%The mathematics behind the specific square root implementation used here
%are described in [2].
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%March 2015, David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
if(isempty(FJacob))
    FJacob=@(x)numDiff(x,f,xDim);
end

xPred=f(xPrev);
F=FJacob(xPrev);
SPred=tria([F*SPrev,SQ]);
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
