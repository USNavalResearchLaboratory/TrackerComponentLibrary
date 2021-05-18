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
%INPUTS: z The zDimX1 vector measurement.
%       SR The zDimXzDim lower-triangular square root of the measurement
%          covariance matrix in the native coordinate system of the
%          measurement.
%    zPred The zDimXnumComp measurement predictions from the filter.
% otherInfo The intermediate results returned in the otherInfo output of
%          the sqrtDDFMeasPred function.
% innovTrans An optional function handle that computes and optionally
%          transforms the value of the difference between the observation
%          and any predicted points. This is called as innovTrans(a,b) and
%          the default if omitted or an empty matrix is passed is
%          @(a,b)bsxfun(@minus,a,b). This must be able to handle sets of
%          values. For a zDimX1 measurement, either of the inputs could be
%          zDimXN in size while one of the inputs could be zDimX1 in size.
%          This only needs to be supplied when a measurement difference
%          must be restricted to a certain range. For example, the
%          innovation between two angles will be 2*pi if one angle is zero
%          and the other 2*pi, even though they are the same direction. In
%          such an instance, a function handle to the
%          wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%          appropriate parameters should be passed for innovTrans.
%
%OUTPUTS: xUpdate The xDimXnumComp updated state vectors.
%         SUpdate The updated xDimXxDimXnumComp lower-triangular square-
%                 root state covariance matrices.
%      innov, Szz The zDimXnumComp innovations and the zDimXzDimXnumComp
%                 square-root innovation covariance matrix are returned in
%                 case one wishes to analyze the consistency of the
%                 estimator or use those values in gating or likelihood
%                 evaluation.
%               W The xDimXzDimXnumComp gain used in the update.
%
%See the comments to the function sqrtDDFMeasPred for an example of usage
%of this function. See the comments to sqrtDDFUpdate for more information
%on the algorithm.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

algorithm=otherInfo.algorithm;
SPred=otherInfo.SPred;
Pxz=otherInfo.Pxz;
Szx=otherInfo.Szx;
Szx2=otherInfo.Szx2;
xPred=otherInfo.xPred;

xDim=size(SPred,1);
numComp=size(SPred,3);
zDim=size(z,1);

xUpdate=zeros(xDim,numComp);
SUpdate=zeros(xDim,xDim,numComp);
innov=zeros(zDim,numComp);
Szz=zeros(zDim,zDim,numComp);
W=zeros(xDim,zDim,numComp);
for k=1:numComp
    if(algorithm==2)%Second order DDF
        %Equation 71 in [1] (Equation 110 in [2])
        Szz(:,:,k)=tria([Szx(:,:,k),SR,Szx2(:,:,k)]);
    else
        %Equation 61 in [1] (Equation 100 in [2])
        Szz(:,:,k)=tria([Szx(:,:,k),SR]);
    end

    %Equation 64 in [1] (Equation 103 in [2])
    W(:,:,k)=(Pxz(:,:,k)/Szz(:,:,k)')/Szz(:,:,k);

    %The innovation
    innov(:,k)=innovTrans(z,zPred(:,k));

    %Equation 65 in [1] (Equation 104 in [2]).
    xUpdate(:,k)=xPred(:,k)+W(:,:,k)*innov(:,k);

    if(algorithm==2)%Second order DDF
        %Equation 115 in [2].
        SUpdate(:,:,k)=tria([SPred(:,:,k)-W(:,:,k)*Szx(:,:,k),W(:,:,k)*SR,W(:,:,k)*Szx2(:,:,k)]);
    else
        %Equation 106 in [2]
        SUpdate(:,:,k)=tria([SPred(:,:,k)-W(:,:,k)*Szx(:,:,k),W(:,:,k)*SR]);
    end
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
