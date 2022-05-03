function [xPred, SPred]=sqrtDiscKalPred(xPrev,SPrev,F,SQ,u)
%SQRTDISCKALPRED Perform the discrete-time prediction step that comes with 
%                the square root implementation of the standard linear
%                Kalman filter with additive process noise.
%
%INPUTS: xPrev The xDimX1 state estimate at the previous time-step.
%        SPrev The xDimXxDim lower-triangular square root of the state
%              covariance matrix at the previous time-step.
%            F An xDimXxDim state transition matrix.
%           SQ The xDimXxDim lower-triangular square root of the process
%              noise covariance matrix.
%            u An optional xDim X1 vector that is the control input. If
%              omitted, no control input is used.
%
%OUTPUTS: xPred The xDimX1 predicted state estimate.
%         SPred The xDimXxDim lower-triangular square root of the
%               predicted state covariance estimate.
%
%The mathematics behind the specific square root implementation used here
%are described in [1], with a flow chart given in Appendix G.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(u))
        u=0; 
    end

%Perform the prediction.
    xPred=F*xPrev+u;
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
