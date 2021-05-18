function [ySqrtUpdate,PInvSqrtUpdate]=sqrtInfoFilterUpdate(ySqrtPred,PInvSqrtPred,z,SR,H)
%%SQRTINFOFILTERUPDATE Perform the measurement update step in a linear
%                      square-root information filter.
%
%INPUTS: ySqrtPred The xDimX1 predicted square root information state. The
%                  predicted information state is always PInvSqrtPred times
%                  the predicted target state estimate.
%     PInvSqrtPred The predicted inverse square root information matrix.
%                  If P is the covariance matrix of a Gaussian state x,
%                  then P=PSqrt*PSqrt' and PInvSqrtPred=inv(PSqrt). This
%                  can be either upper triangular or lower triangular.
%                z The zDim X 1 vector measurement.
%               SR The zDim X zDim lower-triangular square root of  
%                  the measurement covariance matrix. This matrix must be
%                  invertible.
%                H The zDim X xDim measurement matrix for a linear
%                  measurement model. That is z=H*x+w, where w is
%                  measurement noise having covariance matrix R.
%
%OUTPUTS: ySqrtUpdate The xDimX1 updated information state vector.
%      PInvSqrtUpdate The xDimXxDim upper-triangular updated inverse square
%                     root information matrix.
%
%The algorithm is that of the linear square root information filter that
%is described in the paper [1]. However, it has been implemented without
%the additions to allow for noise-free measurement components. The same
%steps are described in the beginning of [2], though this implementation
%considers a decentralized application.
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
%The linear measurement equation addressed by this function is 
%z=H*x_t+v
%where x is the true value of the state, H is the measurement matrix and v
%is Gaussian noise having zero mean and covariance matrix R. For the
%purpose of the square root information filter, instead of R, one uses
%RInvSqrt=inv(chol(R,'lower'))
%
%REFERENCES:
%[1] M. L. Psiaki, "Square-root information filtering and fixed-interval
%    smoothing with singularities," in Proceedings of the American Control
%    Conference, Philadelphia, PA, Jun. 1998, pp. 2744-2748.
%[2] G. J. Bierman and M. R. Belzer, "A decentralized square root
%    information filter/smoother," in Proceedings of the 24th Conference on
%    Decision and Control, Fort Lauderdale, FL, Dec. 1985, pp. 1902-1905.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(ySqrtPred,1);

%A is the matrix on the left-hand side of Equation (12) in Huang's paper.
A=[PInvSqrtPred, ySqrtPred;
   SR\H,   SR\z];

%For an mXn matrix A, we are finding A=QDecomp*T, where T is an mXn
%upper-triangular matrix and QDecomp, which we do not care about, is an mXm
%unitary matrix, meaning QDecomp'*QDecomp==QDecomp*QDecomp'==eye(m,m).
[~,T]=qr(A);
PInvSqrtUpdate=T(1:xDim,1:xDim);
ySqrtUpdate=T(1:xDim,xDim+1);
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
