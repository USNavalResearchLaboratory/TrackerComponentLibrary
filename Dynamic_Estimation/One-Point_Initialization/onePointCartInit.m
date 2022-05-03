function [x,P]=onePointCartInit(zCart,SRCart,higherDerivStdDev,matType)
%%ONEPOINTCARTINIT This function implements single-point initialization for
%              target states that consist of components of position,
%              velocity, acceleration, etc. This function initializes
%              tracks from a single measurement by setting the position
%              components to Cartesian measurement value with its
%              associated covariance and then setting the diagonal elements
%              of the rest of the components based on fixed standard
%              deviations. For example, in [1] it is suggested to make the
%              standard deviation for velocity vMax/3 and in Chapter 3.2.2
%              of [2], vMax/2. Similar ad-hoc values could be used for
%              higher moments. This function does not use Doppler/ range
%              rate.
%
%INPUTS: zCart A zDimXnumMeas set of Cartesian measurements for which
%              single-point differencing should be used to start tracks.
%       SRCart If matType is omitted or is 0, then this is a
%              zDimXzDimXnumMeas set of lower-triangular square root
%              covariance matrices associated with the measurements in
%              zCart. If all of the matrices are the same, then a single
%              zDimXzDim matrix can be passed. If matType=1, then this is a
%              set of covariance matrices.
% higherDerivStdDev A numMomentsX1 or 1XnumMoments vector containing the
%              standard deviations to use for each of the moments
%              (position, velocity, etc) that cannot be estimated from the
%              data. As mentioned in [1] and in and in Chapter 3.2.2
%              of [2], for velocity, this might be vMax/sqrt(2) or
%              vMax/sqrt(3).
%      matType An optional input specifying whether SRCart is a set of
%              lower-triangular square roots of the covariance matrix, or
%              whether it is the set of covariance matrices. Possible
%              values are:
%              0 (The default if omitted or an empty matrix is passed)
%                SRCart holds lower-triangular square root covariance
%                matrices.
%              1 SRCart holds covariance matrices.
%
%OUTPUTS: x The xDimXnumMeas set of target state estimates. All
%           non-position components are zero. xDim=zDim*(numMoments+1). The
%           components are arranged position, velocity, acceleration, etc.
%           For example, [x;y;z;xDot;yDot;zDot].
%         P The xDimXxDimXnumMeas set of initial target state covariance
%           matrices associated with x.
%
%One-point differencing is discussed in [1] and Chapter 3.2.2 of [2].
%
%REFERENCES:
%[1] M. Mallick and B. La Scala, "Comparison of single-point and two-point
%    difference track initiation algorithms using position measurements". 
%    Acta Automatica Sinica, vol.34, no. 3, pp 258-265, Mar. 2008.
%[2] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(matType))
    matType=0; 
end

zDim=size(zCart,1);
numMeas=size(zCart,2);

if(size(SRCart,3)==1)
    SRCart=repmat(SRCart,1,1,numMeas);
end

numMoments=length(higherDerivStdDev);
xDim=zDim*(numMoments+1);

x=zeros(xDim,numMeas);
P=zeros(xDim,xDim,numMeas);

x(1:zDim,:)=zCart;
switch(matType)
    case 0
        for curMeas=1:numMeas
            P(1:zDim,1:zDim,curMeas)=SRCart(:,:,curMeas)*SRCart(:,:,curMeas)';
        end
    case 1
        for curMeas=1:numMeas
            P(1:zDim,1:zDim,curMeas)=SRCart(:,:,curMeas);
        end
    otherwise
        error('Unknown matrix type specified.')
end

sel=(zDim+1):xDim;
P(sel,sel,:)=repmat(kron(diag(higherDerivStdDev(:).^2),eye(zDim,zDim)),1,1,numMeas);
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
