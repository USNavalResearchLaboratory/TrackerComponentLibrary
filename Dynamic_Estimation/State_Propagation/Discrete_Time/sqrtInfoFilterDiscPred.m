function [ySqrtPred, PInvSqrtPred,RwPred,RwxPred]=sqrtInfoFilterDiscPred(ySqrtPrev,PInvSqrtPrev,F,SQ,u,Gamma)
%%SQRTINFOFILTERDISCPRED Perform the discrete-time prediction step that  
%                   comes with the standard linear square root information
%                   filter with additive process noise.
%
%INPUTS: ySqrtPrev The xDimX1 square root information state that is to be
%                  propagated. The previous information state is always
%                  PInvSqrtPrev times the previous target state estimate.
%     PInvSqrtPrev The previous inverse square root information matrix.
%                  If P is the covariance matrix of a Gaussian state x,
%                  then P=PSqrt*PSqrt' and PInvSqrtPrev=inv(PSqrt). This
%                  can be either upper triangular or lower triangular.
%                F An invertible xDim X xDim state transition matrix.
%               SQ The lower-triangular square root of the process noise
%                  covariance matrix. In the standard linear dynamic
%                  model, one has
%                  x(k)=F(k-1)*x(k-1)+u(k-1)+noise where x is the state
%                  and Q is the covariance matrix of the noise. However,
%                  if the covariance matrix of the noise in such equation
%                  is singular, then one should rewrite it as
%                  x(k)=F(k-1)*x(k-1)+Gamma*(u(k-1)+noise)
%                  where Gamma is some matrix and the covariance matrix of
%                  the untransformed noise is not singular. This SQ is the
%                  lower triangular square root of the noise in the linear
%                  dynamic equation.
%                u An optional SQDim X1 vector that is the control input.
%                  If omitted or an empty matrix is passed, a zero control
%                  input (no control input) is used. Note that this control
%                  input differs from formulations in other filters in that
%                  it is affected by the matrix Gamma, if present, which
%                  also transforms the process noise.
%            Gamma An optional matrix that transformes the process noise
%                  and the control input to the state domain if the process
%                  noise covariance matrix is singular, as discussed for
%                  the input SQ. If this is omitted or an empty matrix is
%                  passed, then an identity matrix is used (i.e. there is
%                  no Gamma).
%
%OUTPUTS: ySqrtPred The xDim X 1 predicted square root information state
%                   vector.
%      PInvSqrtPred The predicted xDim X xDim inverse square root state
%                   covariance matrix, which is upper-triangular.
%    RwPred,RwxPred Noise matrices which can be saved for use in smoothing
%
%Given a Gaussian predicted state with mean x and covariance matrix P, the
%information state is
%y=inv(P)*x
%and the information matrix is
%PInv=inv(P)
%On the other hand, the square root information state is 
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
%Note that if there is no process noise or control input, then the
%propagation can be reduced to
%ySqrtPred=ySqrtPrev
%PInvSqrtPred=PInvSqrtPrev/F;
%
%The square root information filter is discussed in Chapter 5, 9 and 10 of
%[1].
%
%The derivation is as follows:
%With a state propagation step of
%x_{k+1}=F*x_{k}+Gamma*(u+w)
%where w is process noise. This can be rewritten as
%x_{k}=F^{-1}*(x_{k+1}-Gamma*(u+w))
%We define
%P=S*S'%Prior state covariance factorization
%Q=W*W'%Process noise covariance factorization
%and the whitened total input is uTotal=inv(W)*u+inv(W)*v;
%The we can write a cost function of
%J=norm(S^{-1}*(F^{-1}*(x_{k}-Gamma*u)-ySqrtPrev)^2+norm(W^{-1}*u-uTotal)^2
%where the input u is known and the terms F^{-1}*x_{k}=ySqrtPred and uTotal
%are unknown. After using a Gram-Schmidt-decomposition to diagonalize the
%cost function in matrix form, one gets a matrix form that can be solved
%for the unknowns. That is the form used for the solution in this function.
%
%REFERENCES:
%[1] G. J. Bierman, Factorization Methods for Discrete Sequential
%    Estimation. Mineola, NY: Dover Publications, Inc., 1977.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(ySqrtPrev,1);
    QDim=size(SQ,1);

    if(nargin<6||isempty(Gamma))
        Gamma=eye(xDim);
    end
    
    if(nargin<5||isempty(u))
       u=zeros(QDim,1); 
    end
    PInvSqrtTilde=PInvSqrtPrev/F;
    A=[inv(SQ),                zeros(QDim,xDim), SQ\u;
       -PInvSqrtTilde*Gamma,   PInvSqrtTilde,   ySqrtPrev];
   
   [~,T] = qr(A);

   ySqrtPred=T((QDim+1):end,end);
   PInvSqrtPred=T((QDim+1):end,(end-xDim):(end-1));
   if(nargout>2)
       RwPred=T(1:QDim,1:QDim);
       RwxPred=T(1:QDim,(QDim+1):(end-1));
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
