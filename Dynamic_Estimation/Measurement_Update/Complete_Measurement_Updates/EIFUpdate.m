function [yUpdate,PInvUpdate,innov]=EIFUpdate(yPred,PInvPred,z,RInv,h,HJacob,innovTrans)
%%EIFUPDATE Perform the measurement update step in the first-order Extended
%           Information Filter (EIF).
%
%INPUTS: yPred The xDimX1 predicted information state. The information
%              state is the inverse covariance matrix times the target
%              state.
%     PInvPred The xDimXxDim inverse of the predicted state covariance
%              matrix.
%            z The zDim X 1 vector measurement.
%         RInv The zDim X zDim inverse of the measurement covariance matrix
%              in the native coordinate system of the measurement.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%       HJacob A function handle for the measurement Jacobian matrix that
%              takes the target state as a parameter or the zDim X xDim
%              measurement Jacobian matrix itself. Each column is the
%              derivative vector of h with respect to the corresponding
%              element of x. If not supplied or an empty matrix is passed,
%              then HJacob will be found using numerical differentiation
%              via the numDiff function with default parameters.
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
%
%OUTPUTS: yUpdate The xDim X 1 updated (posterior) information state
%                 vector.
%      PInvUpdate The updated xDim X xDim inverse state covariance matrix.
%           innov The zDimX1 innovation. This is the difference between the
%                 measurement and the predicted measurement. This is
%                 sometimes used for gating measurements.
%
%This is an implementation of the measurement update step of an extended
%information filter as described in [1].
%
%The optional parameter innovTrans is not described in
%[1], but allow for possible modifications to the filter as
%described in [2]. The parameter has been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] K. P. B. Chandra, D.-W. Gu, and I. Postlethwaite, "Square root
%    cubature information filter," IEEE Sensors Journal, vol. 13, no. 2,
%    pp. 750-758, Feb. 2013.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

if(nargin<7||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

%Extract the state
xPred=lsqminnorm(PInvPred,yPred);

if(isa(HJacob,'function_handle'))
    H=HJacob(xPred);
else
    H=HJacob;%If the Jacobian matrix was directly given.
end

zPred=h(xPred);
innov=innovTrans(z,zPred);

%Equation 8
i=H'*RInv*(innov+H*xPred);

%Equation 9
I=H'*RInv*H;

%Equation 6
yUpdate=yPred+i;

%Equation 7
PInvUpdate=PInvPred+I;

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
