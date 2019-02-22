function [ySqrtPred,PInvSqrtPred,RwPred,RwxPred]=ESRIFDiscPred(ySqrtPrev,PInvSqrtPrev,f,FJacob,SQ,u,Gamma)
%%ESRIFDISCPRED Perform the discrete-time prediction step that comes with
%               the extended square root information filter with additive
%               process noise.
%
%INPUTS: ySqrtPrev The xDimX1 square root information state that is to be
%                  propagated. The previous information state is always
%                  PInvSqrtPrev times the previous target state estimate.
%     PInvSqrtPrev The previous inverse square root information matrix.
%                  If P is the covariance matrix of a Gaussian state x,
%                  then P=PSqrt*PSqrt' and PInvSqrtPrev=inv(PSqrt). This
%                  can be either upper triangular or lower triangular.
%                f A function handle for the state transition function
%                  that takes the state as its parameter.
%           FJacob A function handle for calculating the xDim X xDim state
%                  transition matrix. If an empty matrix is passed, then
%                  FJacob will be found using numerical differentiation 
%                  via the numDiff function with default parameters.
%               SQ The lower-triangular square root of the process noise
%                  covariance matrix. In the standard linear dynamic
%                  model, one has
%                  x(k)=F(k-1)*x(k-1)+u(k-1)+noise where x is the state
%                  and Q is the covariance matrix of the noise. However,
%                  if the covariance matrix of the noise in such equation
%                  is singular, then one should rewrite it as
%                  x(k)=F(k-1)*x(k-1)+u(k-1)+Gamma*noise
%                  where gamma is some matrix and the covariance matrix of
%                  the untransformed noise is not singular. This SQ is the
%                  lower triangular square root of the noise in the linear
%                  dynamic equation.
%                u An optional xDim X1 vector that is the control input.
%                  If omitted, a zero control input (no control input) is
%                  used.
%            Gamma An optional matrix that transforms the process noise
%                  to the state domain if the process noise covariance
%                  matrix is singular, as disccused for the input SQ. If
%                  this is omitted an identity matrix is used (i.e. there
%                  is no Gamma).
%
%OUTPUTS:ySqrtPred The xDim X 1 predicted square root information state
%                  vector.
%     PInvSqrtPred The predicted xDim X xDim inverse square root state
%                  covariance matrix, which is upper-triangular.
%   RwPred,RwxPred  Noise matrices which can be saved for use in smoothing
%
%The time prediction step is taken from the algorithmic implementation
%of [1] described in chapters V and VI.
%
%Given a Gaussian predicted state with mean x and covariance matrix P, the
%square root information state is
%ySqrt=PInvSqrt*x
%where
%PSqrt=inv(PInvSqrt)
%and
%P=PSqrt*PSqrt';
%The matrix PInvSqrtPred can be upper or lower triangular, when supplied to
%this function. For example, a lower-triangular matrix can be obtained
%using PInvSqrt=inv(chol(P,'lower')). However, the output of this function,
%PInvSqrtUpdate is always upper triangular.
%
%REFERENCES:
%[1] G. J. Bierman, "Factorization Methods for Discrete Sequential
%    Estimation. Academic Press, New York, 1977.
%
%February 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(ySqrtPrev,1);
if(nargin<6||isempty(u))
    u=zeros(xDim,1);    
end
if(nargin<7||isempty(Gamma))
    Gamma=eye(xDim);
end
if(isempty(FJacob))
    FJacob=@(x)numDiff(x,f,xDim);
end

xPrev=PInvSqrtPrev\ySqrtPrev;

F=FJacob(xPrev);
PInvSqrtTilde=PInvSqrtPrev/F;

%The mapping matrix found on p.121 of Bierman
A=[inv(SQ),                zeros(xDim,xDim), SQ\u;
    -PInvSqrtTilde*Gamma,   PInvSqrtTilde,   ySqrtPrev];
[~,T] = qr(A);

PInvSqrtPred=T((xDim+1):end,(end-xDim):(end-1));

%Since f may be nonlinear, the information state output described in
%Bierman is not valid, and the predicted state should be calculated
%independently.
ySqrtPred=PInvSqrtPred*f(xPrev);

if(nargout>2)
    RwPred=T(1:xDim,1:xDim);
    RwxPred=T(1:xDim,xDim+1:end-1);
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
