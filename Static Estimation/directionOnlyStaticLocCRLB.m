function CRLB=directionOnlyStaticLocCRLB(t,lRx,M,SRLocal,xi,w)
%%DIRECTIONONLYSTATICLOCCRLB Obtain the Cramér-Rao lower bound (CRLB) for
%             localizing a target given measurements in the form of
%             direction cosines (two components of a unit vector in 3D, or
%             one component in 2D) in the local coordinates systems of the
%             receivers. The target is assumed to be in front of the
%             receivers, so the final component is uniquely defined to
%             obtain a unit vector and is positive. The measurements are
%             approximated as being corrupted with Gaussian noise having a
%             known covariance matrix in the local coordinate systems. The
%             assumption is that the measurements are sufficiently far from
%             the edge of the observable region/ the noise is sufficiently
%             low that the measured components will never exceed valid
%             values to be one (in 2D) or two (in 3D) components of unit
%             vectors.
%
%INPUTS:  t The 2X1 (in 2D) or 3X1 (in 3D) location of the target in global
%           coordinates. The target must be "in front" of all of the
%           sensors. That means that when a unit vector from the receiver
%           to the target is transformed into the local coordinate system
%           of each receiver, the third coordinate (the z-coordinate, which
%           is the boresight direction) must be positive.
%       lRx The 2XnumMeas or 3XnumMeas set of Cartesian locations of the
%           sensors producing the direction vector measurements.
%         M The 2X2XnumMeas or 3X3XnumMeas set of rotation matrices to
%           rotate from global coordinates into the local cordinate system
%           of each sensor. Note that the boresight (pointing direction) of
%           the sensor is the local z-axis. The u-v components in the
%           measurements are along local x and y axes.
%   SRLocal A 2X2XnumMeas (for 3D) or 1X1XnumMeas (for 2D) set of
%           lower-triangular square roots of the covariance matrices for
%           the local measurement components. These must be invertible.
%      xi,w Optionally, cubature points and weights for integration over a
%           multivariaute Gaussian PDF can be supplied. The points must be
%           (numDim-2)XnumMeas-dimensional. If these are omitted or empty
%           matrices are passed, then the default set of fifth-order
%           cubature points using the function fifthOrderCubPoints are
%           used. The points are necessary for evaluating the expected
%           value needed to compute the Fisher information matrix.
%
%OUTPUTS: CRLB A 2X2 or 3X3 approximation of the CRLB. This is a
%              lower-bound on the covariance matrix of an unbiased
%              estimator of the target location given the measurements. The
%              CRLB is the inverse of the Fisher information matrix.
%
%The CRLB algorithm is implemented based on [1]. The CRLB is also discussed
%in Section 2.7.2 of [2].
%
%EXAMPLE:
%4 sensors on Hawaiian islands
% lRxLatLonAlt=[[20.72529087;-156.1164093;1618],...
%               [20.74070316;-155.98731995;36],...
%               [20.25189031;-155.80604553;53],...
%               [20.1152606;-155.54855347;187]];
% lRxLatLonAlt(1:2,:)=lRxLatLonAlt(1:2,:)*(pi/180);%Convert to radians.
% lRx=ellips2Cart(lRxLatLonAlt);%Convert to Cartesian
% targetLatLonAlt=[20.62250226;-155.42495728;8000];%Middle scenario
% targetLatLonAlt(1:2)=targetLatLonAlt(1:2)*(pi/180);%Convert to radians.
% t=ellips2Cart(targetLatLonAlt);%Convert to Cartesian
% 
% %Two radars pointed East, Two radars pointed North.
% az12=90*(pi/180);%East
% el=15*(pi/180);%Radar elevation, radians.
% %M will be rotation matrices to go from global to local coordinates.
% M=zeros(3,3,4);
% M(:,:,1)=findRFTransParam(lRxLatLonAlt(:,1),az12,el);
% M(:,:,2)=findRFTransParam(lRxLatLonAlt(:,2),az12,el);
% %The final two radars are pointing North
% az34=0*(pi/180);%North
% M(:,:,3)=findRFTransParam(lRxLatLonAlt(:,3),az34,el);
% M(:,:,4)=findRFTransParam(lRxLatLonAlt(:,4),az34,el);
% 
% %All sensors have the same accuracy in the local coordinate system. About
% %1mrad, when considering points near the boresight.
% SRLocal=repmat(diag([1e-3;1e-3]),1,1,4);
% CRLB=directionOnlyStaticLocCRLB(t,lRx,M,SRLocal)
% %The minimum root-mean-squared error of an estimator is the square root
% %of the trace of the CRLB. 
% minRMSE=sqrt(trace(CRLB))
% %The minRMSE should be about 52.4562.
%
%REFERENCES:
%[1] D. F. Crouse, "Bearings-Only Localization Using Direction Cosines," in
%    Proceedings of the 19th International Conference on Information
%    Fusion, Heidelberg, Germany, July 2016.
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

numDim=size(lRx,1);
numDir=size(lRx,2);

%The total number of variables over which cubature integration will be
%performed. The approximation used is that the u-v measurement components
%(or just u component for 2D estimation) in the local coordinate systems
%are corrupted with Gaussian noise. Thus, there is one fewer dimension than
%the length of the unit vector.
numVar=(numDim-1)*numDir;

if(nargin<5||isempty(xi))%Use default cubature points and weights.
    [xi,w]=fifthOrderCubPoints(numVar);
elseif(size(xi,1)~=numVar)
    error('The cubature points passed have the wrong dimensionality')
end

%u is a set of the true direction vectors in the global coordinate system.
u=bsxfun(@minus,t,lRx);
u=bsxfun(@rdivide,u,sqrt(sum(u.*u,1)));
%Transform the global direction vectors into direction vectors in the local
%coordinate systems:
uRx=zeros(numDim,numDir);
for i=1:numDir
    %M to go from global to local
    uRx(:,i)=M(:,:,i)*u(:,i);
end

if(any(uRx(end,:)<0))
   error('The target should be in front of all of the sensors in their local coordinate systems.')
end

numCubPoints=length(w);

%The covariance matrix to use for the cubature points for all of the
%values.
SRLarge=blkDiagRep(SRLocal);

measVals=uRx(1:(numDim-1),:);

%Transform the cubature points.
xi=bsxfun(@plus,measVals(:),SRLarge*xi);
xi=reshape(xi,numDim-1,numDir*numCubPoints);

%The third dimension in the local coordinate systems must be added to the
%vectors. As the target is in front of all of the sensors, the sign of the
%third coordinate is always positive.
%The abs deals with points near the edge that get pushed into an invalid
%region.
xi=[xi;sqrt(abs(1-sum(xi.*xi)))];
xi=reshape(xi,numDim,numDir,numCubPoints);

%Rotate into global coordinates
for curPoint=1:numCubPoints
    for curDir=1:numDir
        xi(:,curDir,curPoint)=M(:,:,curDir)'*xi(:,curDir,curPoint);
    end
end

%Get the necessary 3D inverse covariance matrices.
RInvGlobal=zeros(3,3,4);%Allocate space
for i=1:numDir
    %The third coordinate provides no additional information locally, so
    %the inverse covariance values for it are set to zero (a singular
    %inverse covariance matrix).
    RInvGlobal(1:(numDim-1),1:(numDim-1),i)=inv(SRLocal(:,:,i)*SRLocal(:,:,i)');
    RInvGlobal(:,:,i)=M(:,:,i)*RInvGlobal(:,:,i)*M(:,:,i)';
end

%Allocate space for the Fisher information matrix.
J=zeros(numDim,numDim);
%Take the expected value voer all of the cubature points.
for curPoint=1:numCubPoints
    J=J+w(curPoint)*logLikelihoodOuterProduct(t,xi(:,:,curPoint),lRx,RInvGlobal);
end

CRLB=inv(J);
end

function J=logLikelihoodOuterProduct(t,u,lRx,RInv)

numDim=size(u,1);
numDir=size(u,2);
J=zeros(numDim,numDim);

for i=1:numDir
    for j=i:numDir
        diffi=t-lRx(:,i);
        diffiMag=norm(diffi);
        diffj=t-lRx(:,j);
        diffjMag=norm(diffj);
        
        Ai=getA(diffi);
        Aj=getA(diffj);
        
        term1=(1/(diffiMag^4*diffjMag^4))*Ai*RInv(:,:,i)*diffi*diffj'*RInv(:,:,j)*Aj;
        term2=(1/(diffiMag^3*diffjMag^3))*Ai*RInv(:,:,i)*u(:,i)*u(:,j)'*RInv(:,:,j)*Aj;
        term3=-(1/(diffiMag^4*diffjMag^3))*Ai*RInv(:,:,i)*diffi*u(:,j)'*RInv(:,:,j)*Aj;
        term4=-(1/(diffiMag^3*diffjMag^4))*Ai*RInv(:,:,i)*u(:,i)*diffj'*RInv(:,:,j)*Aj;
        
        curTerm=2*(term1+term2+term3+term4);
        
        if(i==j)
            J=J+curTerm;
        else
            J=J+curTerm+curTerm';
        end
    end
end

end

function A=getA(tlDiff)
%Get the matrices of partial derivatives.
    numDim=size(tlDiff,1);
    if(numDim==2)
        A=[tlDiff(2)^2,            -tlDiff(1)*tlDiff(2);
           -tlDiff(1)*tlDiff(2),     tlDiff(1)^2];
    else
        A=[tlDiff(2)^2+tlDiff(3)^2, -tlDiff(1)*tlDiff(2),   -tlDiff(1)*tlDiff(3);
           -tlDiff(1)*tlDiff(2),    tlDiff(1)^2+tlDiff(3)^2,-tlDiff(2)*tlDiff(3);
           -tlDiff(1)*tlDiff(3),    -tlDiff(2)*tlDiff(3),   tlDiff(1)^2+tlDiff(2)^2];
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
