function [yPred,PInvPred]=discCubInfoPred(yPrev,PInvPrev,f,Q,xi,w,stateDiffTrans,stateAvgFun,stateTrans)
%%DISCCUBINFOPRED Perform the discrete-time prediction step that comes  
%                 the cubature information filter with additive process
%                 noise. As the cubature information filter as derived in
%                 [1] has no real propagation step, this just extracts the
%                 state, calls discCubKalPred, and converts the result back
%                 to an information state. This function is useful if one
%                 is using an information state because of how the
%                 measurement update in an information filter is performed.
%
%INPUTS: yPrev The xDimX1 information state at the previous time-step. The
%              information state is the inverse covariance matrix times the
%              target state.
%     PInvPrev The xDimXxDim inverse of the state covariance matrix at the
%              previous time-step.
%            f A function handle for the state transition function that
%              takes the state as its parameter.
%            Q The xDimX xDim process noise covariance matrix.
%           xi A (xDim+cDim) X numCubPoints matrix of cubature points. If
%              this and the next parameter are omitted or empty matrices
%              are passed, then fifthOrderCubPoints(xDim+cDim) is used. It
%              is suggested that xi and w be provided to avoid needless
%              recomputation of the cubature points.     
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points.
% stateDiffTrans An optional function handle that takes an xDimXN matrix of
%                N differences between states estimates and transforms them
%               however might be necessary. For example, a state continaing
%               angular components will generally need differences between
%               angular components wrapped to the range +/-pi.
%   stateAvgFun An optional function that given an xDimXN matrix of N state
%               estimates and an NX1 vector of weights, provides the
%               weighted average of the state estimates. This is necessary
%               if, for example, states with angular components are
%               averaged.
%    stateTrans An optional function that takes a state estimate and
%               transforms it. This is useful if one wishes the elements of
%               the state to be bound to a certain domain. For example, if
%               an element of the state is an angle, one should generally
%               want to bind it to the region +/-pi. This is not applied to
%               the output of f.
%
%OUTPUTS: yPred  The xDim X 1 predicted information state vector.
%       PInvPred The predicted xDim X xDim inverse state covariance matrix.
%
%The idea of a cubature information filter is given in Algorithm 1 in [1],
%but the prediction step just converts back to a non-information state and
%performs the update the traditional way before converting back. In other
%words, there is no prediction step specially tailored to an information
%state in that filter.
%
%REFERENCES:
%[1] K. P. B. Chandra, D.-W. Gu, and I. Postlethwaite, "Square root
%    cubature information filter," IEEE Sensors Journal, vol. 13, no. 2,
%    pp. 750-758, Feb. 2013.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
   xi=[]; 
end

if(nargin<5)
   w=[]; 
end

if(nargin<7||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

if(nargin<8||isempty(stateAvgFun))
    stateAvgFun=@(x,w)calcMixtureMoments(x,w);
end

if(nargin<9||isempty(stateTrans))
    stateTrans=@(x)x;
end

PPrev=pinv(PInvPrev);
xPrev=PPrev*yPrev;

[xPred, PPred]=discCubKalPred(xPrev,PPrev,f,Q,xi,w,stateDiffTrans,stateAvgFun,stateTrans);

PInvPred=pinv(PPred);
yPred=PInvPred*xPred;

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
