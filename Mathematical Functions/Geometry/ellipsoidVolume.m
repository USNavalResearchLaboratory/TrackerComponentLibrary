function V=ellipsoidVolume(A,gammaVal,useAInv)
%%ELLIPSOIDVOLUME Return the volume of an ellipsoid, where a point on the
%              ellipsoid satisfies the equation (zp-z)'*A*(zp-z)=gammaVal.
%              This type of volume tends to arise in target tracking when
%              performing gating. In such an instance, inv(A) is a Gaussian
%              covariance matrix associated with the estimate and the
%              ellipsoid is a probability region, where gammaVal determines
%              what amount of probability is in the region. 
%
%INPUTS: A A numDimXnumDimXN set of N positive definite matrices that
%          specify the size and shape of the ellipses or ellipsoids, where
%          a point zp is on the ith ellipse/ ellipsoid if
%          (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gammaVal. Alternatively, if
%          the input useAInv values are passed, it is assumed that this
%          parameter is actually the inverse of the desired matrix, so the
%          region is (zp-z(:,i))'*inv(A(:,:,i))*(zp-z(:,i))=gammaVal
% gammaVal A parameter specifying the size of the ellipse/ ellipsoid.
%          gammaVal must be positive. To specify a probability region of
%          probReg as is commonly used in tracking where A is a Gaussian
%          covariance matrix, one can get gammaVal from
%          ChiSquareD.invCDF(probReg,size(A,1)) The default if this
%          parameter is omitted or an empty matrix is passed is 1.
%  useAInv This parameter indicates whether A should be inverted, as
%          described above. The default if this parameter is omitted or an
%          empty matrix is passed is false.
%
%OUTPUTS: V The NX1 set of scalar volumes of the gating regions specified
%           by the A and gammaVal values.
%
%The formula for the volume of an ellipsoidal gating region is from Chapter
%2.3.2 of [1].
%
%REFERENCES:
%[1] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(useAInv))
    useAInv=false;
end

if(nargin<2||isempty(gammaVal))
    gammaVal=1; 
end

numOut=size(A,3);

n=size(A,1);
cn=(pi^(n/2)/gamma(n/2+1));

V=zeros(numOut,1);
if(useAInv==false)
    for curPoint=1:numOut
        V(curPoint)=cn*sqrt(1/det(1/gammaVal*A(:,:,curPoint)));
    end
else
    for curPoint=1:numOut
        V(curPoint)=cn*sqrt(det(gammaVal*A(:,:,curPoint)));
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
