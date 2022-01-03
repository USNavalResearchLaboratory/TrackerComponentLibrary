function [xUpdate,SUpdate,innov,Szz,W]=sqrtEKFUpdate(xPred,SPred,z,SR,h,HJacob,innovTrans,stateTrans)
%SQRTEKFUPDATE Perform the measurement update step in a square-root version
%              of the first-order Extended Kalman Filter (EKF).
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular square root predicted state
%              covariance matrix.
%            z The zDim X 1 measurement vector.
%           SR The zDim X zDim lower-triangular square root of the
%              measurement covariance matrix in the native coordinate
%              system of the measurement.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%       HJacob A function handle for the measurement Jacobian matrix that
%              takes the target state as a parameter. If not supplied or an
%              empty matrix is passed, then HJacob will be found using
%              numerical differentiation via the numDiff function with
%              default parameters.
%   innovTrans An optional function handle that computes and optionally
%              transforms the value of the difference between the
%              observation and any predicted points. This is called as
%              innovTrans(a,b) and the default if omitted or an empty
%              matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%              able to handle sets of values. For a zDimX1 measurement,
%              either of the inputs could be zDimXN in size while one of
%              the inputs could be zDimX1 in size.  This only needs to be
%              supplied when a measurement difference must be restricted
%              to a certain range. For example, the innovation between two
%              angles will be 2*pi if one angle is zero and the other
%              2*pi, even though they are the same direction. In such an
%              instance, a function handle to the
%              wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%              appropriate parameters should be passed for innovTrans.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one should generally
%              want to bind it to the region +/-pi.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         SUpdate The updated xDim X xDim lower-triangular square root
%                 state covariance matrix.
%      innov, Szz The zDimX1 innovation and the zDimXzDim square root
%                 innovation covariance matrix are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%The first-order EKF is summarized in Figure 10.3.3-1 in Chapter 10.3.3 of 
%[1]. The Joseph-form covariance update given in Chapter 5 of the same book
%is used for improved numerical stability. The mathematics behind the
%specific square root implementation used here are those used in the
%standard Kalman filter, which are described in [2] with a flow chart given
%in Appendix G.
%
%The optional parameter innovTrans is not described in the above
%reference, but allows for possible modifications to the filter as
%described in [3].
%
%The parameters have been added to allow the filter to be used with
%angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[3] David F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering
%    with angular measurement models," in Proceedings of the 18th
%    International Conference on Information Fusion, Washington, D.C.,
%    6-9 Jul. 2015.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%Updated with additional inputs, August 2015 David F. Crouse,
%                               Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zDim=size(z,1);

if(nargin<8||isempty(stateTrans))
    stateTrans=@(x)x;
end

if(nargin<7||isempty(innovTrans))
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

if(nargin<6||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

zPred=h(xPred);

H=HJacob(xPred);
Pxz=SPred*SPred'*H';

Szz=tria([H*SPred,SR]);
W=(Pxz/Szz')/Szz;

innov=innovTrans(z,zPred);
xUpdate=stateTrans(xPred+W*innov);

temp=W*H;
SUpdate=tria([(eye(size(temp))-temp)*SPred,W*SR]);

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
