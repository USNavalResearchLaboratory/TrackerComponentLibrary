function [xUpdate,SUpdate,innov,Szz,W]=sqrtDDFUpdateWithPred(z,SR,zPred,otherInfo,innovTrans)
%%SQRTDDFUPDATEWITHPRED Given the output of the measurement prediction step
%           from sqrtDDFMeasPred and a measurement, complete the
%           measurement update step of the square-root version of the
%           central difference filter (CDF), or the first or second order
%           divided difference filter (DDF) with additive measurement
%           noise. Separating the measurement prediction step from the rest
%           of the update step can make the creation of multiple
%           measurement association hypotheses from a single target
%           prediction more efficient. The full measurement update function
%           is sqrtDDFUpdate.
%
%INPUTS: z The zDim X 1 vector measurement.
%       SR The zDim X zDim lower-triangular square root of the measurement
%          covariance matrix in the native coordinate system of the
%          measurement.
%    zPred The zDimX1 measurement prediction from the filter.
% otherInfo The intermediate results returned in the otherInfo output of
%          the sqrtDDFMeasPred function.
%   innovTrans An optional function handle that transforms the value of the
%          difference between the observation and any predicted points.
%          This must be able to handle sets of differences. For a zDim 
%          measurement, this must be able to handle a zDimXN matrix of N
%          differences. This only needs to be supplied when a measurement
%          difference must be restricted to a certain range. For example,
%          the innovation between two angles will be 2*pi if one angle is
%          zero and the other 2*pi, even though they are the same
%          direction. In such an instance, a function handle to the
%          wrapRange function with the appropriate parameters should be
%          passed for innovTrans.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         SUpdate The updated xDim X xDim lower-triangular square root
%                 state covariance matrix.
%      innov, Szz The zDimX1 innovation and the zDimXzDim square root
%                 innovation covariance matrix are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update.
%
%See the comments to the function sqrtDDFMeasPred for an example of usage
%of this function. See the comments to sqrtDDFUpdate for more information
%on the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(x)x;
end

algorithm=otherInfo.algorithm;
SPred=otherInfo.SPred;
Pxz=otherInfo.Pxz;
Szx=otherInfo.Szx;
Szx2=otherInfo.Szx2;
xPred=otherInfo.xPred;

if(algorithm==2)%Second order DDF
    %Equation 71 in [1] (Equation 110 in [2])
    Szz=tria([Szx,SR,Szx2]);
else
    %Equation 61 in [1] (Equation 100 in [2])
    Szz=tria([Szx,SR]);
end

%Equation 64 in [1] (Equation 103 in [2])
W=(Pxz/Szz')/Szz;

%The innovation
innov=innovTrans(z-zPred);

%Equation 65 in [1] (Equation 104 in [2]).
xUpdate=xPred+W*innov;

if(algorithm==2)%Second order DDF
    %Equation 115 in [2].
    SUpdate=tria([SPred-W*Szx,W*SR,W*Szx2]);
else
    %Equation 106 in [2]
    SUpdate=tria([SPred-W*Szx,W*SR]);
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
