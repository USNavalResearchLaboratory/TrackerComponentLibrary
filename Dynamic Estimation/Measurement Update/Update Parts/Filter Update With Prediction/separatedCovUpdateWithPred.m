function [xUpdate,LUpdate,PUpdate,innov,Pzz,W]=separatedCovUpdateWithPred(z,R,zPred,PzPred,otherInfo,LPred,c)
%%SEPARATEDCOVUPDATEWITHPRED Given the output of the measurement prediction
%           step from separatedCovMeasPred and a measurement, complete the
%           measurement update step of the separated covariance filter.
%           Separating the measurement prediction step from the rest of the
%           update step can make the creation of multiple measurement
%           association hypotheses from a single target prediction more
%           efficient. The full measurement update function is
%           separatedCovUpdate.
%
%INPUTS: z The zDim X 1 vector measurement.
%        R The zDim X zDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimX1 measurement prediction from the filter.
%   PzPred The zDimXzDim covariance matrix associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the separatedCovMeasPred function.
%    LPred The xDimXzDim predicted delay vector (defined before Equation 12
%          in [1]). The use of multiple columns represents a choice in how
%          the algorithm was generalized to multiple dimensions.
%        c The confidence region under consideration by the filter. 0<c<1.
%          If this parameter is omitted or an empty matrix is passed, the
%          default value of c=0.99 is used.
%
%OUTPUTS: xUpdate The xDimX1 updated state vector.
%         LUpdate The xDimXzDim updated delay vector.
%         PUpdate The xDimXxDim covariance matrix of the state estimate.
%                 This is a combination of LUpdate and PUpdate.
%      innov, Pzz The zDimX1 innovation and a zDimXzDim matrix S that is
%                 akin to an innovation covariance matrix are returned in
%                 case one wishes to analyze the consistency of the
%                 estimator or use those values in gating or likelihood
%                 evaluation.
%               W The gain used in the update.
%
%See the comments to the function separatedCovMeasPred for an example of
%usage of this function. See the comments to separatedCovUpdate for more
%information on the algorithm.
%
%REFERENCES:
%[1] G. J. Portmann, J. R. Moore, and W. G. Bath, "Separated covariance
%    filtering," in Proceedings of the IEEE International Radar Conference,
%    Arlington, VA, 7-10 May 1990, pp. 456-460.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(c))
    c=0.99; 
end

H=otherInfo.H;
TPred=otherInfo.TPred;
xPred=otherInfo.xPred;
Pxz=otherInfo.Pxz;

Pzz=PzPred+c^2*R;

W=Pxz/Pzz;%The gain

innov=z-zPred;%The innovation
xUpdate=xPred+W*innov;

xDim=size(xPred,1);
diff=eye(xDim,xDim)-W*H;

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
