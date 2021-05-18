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
%INPUTS: z The zDimX1 vector measurement.
%        R The zDimXzDim measurement covariance matrix in the native
%          coordinate system of the measurement.
%    zPred The zDimXnumComp measurement predictions from the filter.
%   PzPred The zDimXzDimXnumComp covariance matrices associated with zPred.
% otherInfo The intermediate results returned in the otherInfo output of
%          the separatedCovMeasPred function.
%    LPred The xDimXzDimXnumComp predicted delay vectors (defined before
%          Equation 12 in [1]). The use of multiple columns represents a
%          choice in how the algorithm was generalized to multiple
%          dimensions.
%        c The confidence region under consideration by the filter. 0<c<1.
%          If this parameter is omitted or an empty matrix is passed, the
%          default value of c=0.99 is used.
%
%OUTPUTS: xUpdate The xDimXnumComp updated state vectors.
%         LUpdate The xDimXzDimXnumComp updated delay vectors.
%         PUpdate The xDimXxDimXnumComp covariance matrices of the state
%                 estimates. This is a combination of LUpdate and PUpdate.
%      innov, Pzz The zDimXnumComp innovations and a zDimXzDimXnumComp
%                 set of matrices Pzz that are akin to an innovation
%                 covariance matrices are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDimXnumComp gains used in the update.
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

xDim=size(xPred,1);
numComp=size(xPred,2);
zDim=size(z,1);

xUpdate=zeros(xDim,numComp);
LUpdate=zeros(xDim,zDim,numComp);
PUpdate=zeros(xDim,xDim,numComp);
innov=zeros(zDim,numComp);
Pzz=zeros(zDim,zDim,numComp);
W=zeros(xDim,zDim,numComp);
for k=1:numComp
    Pzz(:,:,k)=PzPred(:,:,k)+c^2*R;
    W(:,:,k)=Pxz(:,:,k)/Pzz(:,:,k);%The gain

    innov(:,k)=z-zPred(:,k);%The innovation
    xUpdate(:,k)=xPred(:,k)+W(:,:,k)*innov(:,k);

    xDim=size(xPred,1);
    diff=eye(xDim,xDim)-W(:,:,k)*H;

    LUpdate(:,:,k)=diff*LPred(:,:,k);
    TUpdate=diff*TPred(:,:,k);

    PUpdate(:,:,k)=(TUpdate-LUpdate(:,:,k)*LUpdate(:,:,k)')/c^2;
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
