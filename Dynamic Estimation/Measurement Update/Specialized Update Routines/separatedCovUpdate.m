function [xUpdate,LUpdate,PUpdate]=separatedCovUpdate(xPred,LPred,TPred,z,R,H,c)
%%SEPARATEDCOVUPDATE Perform the measurement update step in the separated
%                    covariance filter. This filter takes a maximum assumed
%                    acceleration (or other moment) for the target and
%                    provides the optimum estimate in terms of a cost
%                    function that trades off between estimation accuracy
%                    and estimator delay. The filter assumes a linear
%                    measurement model.
%
%INPUTS: xPred The xDimX1 predicted target state.
%        LPred The xDimXzDim predicted delay vector (defined before
%              Equation 12 in [1]). The use of multiple columns represents
%              a choice in how the algorithm was generalized to multiple
%              dimensions.
%        TPred The xDimXxDim total error matrix. This is a combination of
%              errors due to measurement noise and filter lag.
%            z The zDim X 1 vector measurement. zDim should be the same as
%              the number of position components in the state.
%            R The zDim X zDim measurement covariance matrix.
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%            c The confidence region under consideration by the filter.
%              0<c<1. If this parameter is omitted or an empty matrix is
%              passed, the default value of c=0.99 is used.
%
%OUTPUTS: xUpdate The xDimX1 updated state vector.
%         LUpdate The xDimXzDim updated delay vector.
%         PUpdate The xDimXxDim covariance matrix of the state estimate.
%                 This is a combination of LUpdate and PUpdate
%
%The equations for the algorithm are given in Table 1 of [1]. The algorithm
%is presented in 1D for a state consisting of position and velocity.
%However, the equations are given in vector form and can thus be used in
%multiple dimensions, which is done here.
%
%The filter in 1 is 1D. To generalize the filter to 3D, L is redefined so
%that each column contains the lag for one particular dimension of motion.
%To keep the solution the same as in Table 1, L becomes a matrix where the
%elements in each column that do not correspond to components for that
%dimensions of motion are zero.
%
%In [1], no clear method of initializing this type of tracking filter is
%provided. A simple way to initialize the filter would be to use two
%Cartesian converted measurements to obtain a state estimate and covariance
%PInit as one would do with a normal Kalman filter (one could, for example,
%use the KalmanFIRSmoother function) and then set TUpdate=c^2*PInit;
%LUpdate=zeros(xDim,zDim); and PUpdate=PInit.
%
%REFERENCES:
%[1] G. J. Portmann, J. R. Moore, and W. G. Bath, "Separated covariance
%    filtering," in Proceedings of the IEEE International Radar Conference,
%    Arlington, VA, 7-10 May 1990, pp. 456-460.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(c))
   c=0.99; 
end

K=TPred*H'/(H*TPred*H'+c^2*R);%The gain

innov=z-H*xPred;%The innovation
xUpdate=xPred+K*innov;

xDim=size(xPred,1);
diff=eye(xDim,xDim)-K*H;

LUpdate=diff*LPred;
TUpdate=diff*TPred;

PUpdate=(TUpdate-LUpdate*LUpdate')/c^2;
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
