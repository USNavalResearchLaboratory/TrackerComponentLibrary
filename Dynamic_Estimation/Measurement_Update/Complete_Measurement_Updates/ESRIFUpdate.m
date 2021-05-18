function [ySqrtUpdate,PInvSqrtUpdate]=ESRIFUpdate(ySqrtPred,PInvSqrtPred,z,SR,h,HJacob,innovTrans)
%%ESRIFUPDATE Perform the measurement update step in a first-order extened
%             square-root information filter.
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
%                h A function handle for the measurement function that
%                  takes the state as its argument.
%           HJacob A function handle for the measurement Jacobian matrix
%                  that takes the target state as a parameter.If an empty
%                  matrix is passed, then HJacob will be found using
%                  numerical differentiation via the numDiff function with
%                  default parameters.
%       innovTrans An optional function handle that computes and optionally
%                  transforms the value of the difference between the
%                  observation and any predicted points. This is called as
%                  innovTrans(a,b) and the default if omitted or an empty
%                  matrix is passed is @(a,b)bsxfun(@minus,a,b). This must
%                  be able to handle sets of values. For a zDimX1
%                  measurement, either of the inputs could be zDimXN in
%                  size while one of the inputs could be zDimX1 in size.
%                  This only needs to be supplied when a measurement
%                  difference must be restricted to a certain range. For
%                  example, the innovation between two angles will be 2*pi
%                  if one angle is zero and the other 2*pi, even though
%                  they are the same direction. In such an instance, a
%                  function handle to the
%                  wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%                  appropriate parameters should be passed for innovTrans.
%
%OUTPUTS: ySqrtUpdate The xDimX1 updated information state vector.
%      PInvSqrtUpdate The xDimXxDim upper-triangular updated inverse square
%                     root information matrix.
%
%The algorithm is that of the extended square root information filter that
%is described in the paper [1], which implements an extended version of the
%square root information filter of [2].
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
%The optional parameter innovTrans is not described in the above
%reference, but allows for possible modifications to the filter as
%described in [3]. The parameters have been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] Psiaki, Mark L., et al. "ALEXIS spacecraft attitude reconstruction
%    with thermal/flexible motions due to launch damage." Journal of
%    Guidance, Control, and Dynamics 20.5 (1997): 1033-1041.
%[2] G. J. Bierman, "Factorization Methods for Discrete Sequential
%    Estimation. Academic Press, New York, 1977.
%[3] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%February 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%Updated with an additional input, August 2015 David F. Crouse,
%                               Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(innovTrans))
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

if(nargin<6||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

xPred=PInvSqrtPred\ySqrtPred;
H=HJacob(xPred);
xDim=size(ySqrtPred,1);

%Adjust measurement to account for nonlinearity in model.
zAdj=H*xPred+innovTrans(z,h(xPred));

%Equation 5 from Psiaki, which uses the adjusted measurement rather than
%the raw measurement as described in Bierman (p.121)
A=[PInvSqrtPred, ySqrtPred;
    SR\H,   SR\zAdj];
[~,T] = qr(A);
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
