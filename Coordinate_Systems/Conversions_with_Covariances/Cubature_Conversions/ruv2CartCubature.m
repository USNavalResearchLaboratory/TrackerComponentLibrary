function [zCart,RCart]=ruv2CartCubature(z,SR,useHalfRange,zTx,zRx,M,xi,w)
%RUV2CARTCUBATURE Use cubature integration to approximate the mean and
%                 covariance matrix of noisy measurements converted from
%                 bistatic r-u-v coordinates into Cartesian coordinates.
%                 For a two-way monostatic conversion, set zTx=zRx to make
%                 the transmitter and receiver collocated. In the
%                 monostatic case, this function is more accurate than 
%                 monostatRuv2CartTaylor, but it is also slower.
%
%INPUTS: z A 3XnumMeas matrix of numMeas vectors to convert. Each has
%          elements [r;u;v], where r is the bistatic range from the
%          transmitter to the target to the receiver, and u and v are
%          direction cosines.
%       SR The 3X3XnumMeas lower-triangular square roots of the measurement
%          covariance matrices for the measurements. If the matrices are
%          all the same, then a single 3X3 matrix can be passed.
% useHalfRange A boolean value specifying whether the bistatic range value
%          has been divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided is false.
%      zTx The 3X1 [x;y;z] location vector of the transmitter in global
%          Cartesian coordinates.  If this parameter is omitted or an empty
%          matrix is passed, then the transmitter is placed at the origin.
%      zRx The 3X1 [x;y;z] location vector of the receiver in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to the local alignment of the receiver. The z
%          axis of the local coordinate system of the receiver is the
%          pointing direction of the receiver. If this matrix is omitted or
%          an empty matrix is passed, then the identity matrix is used.
%       xi A 3XnumCubaturePoints matrix of cubature points for the numeric
%          integration. If this and the final parameter are omitted or
%          empty matrices are passed, then fifthOrderCubPoints is used to
%          generate cubature points.
%        w A numCubaturePointsX1 vector of the weights associated with the
%          cubature points.
%
%OUTPUTS: zCart The approximate means of the PDFs of the Cartesian 
%               converted measurements in [x;y;z] Cartesian coordinates.
%               This is a 3XnumMeas matrix.
%         RCart The approximate 3X3XnumMeas set of covariance matrices of
%               the PDFs of the Cartesian converted measurements.
%
%Details of the conversion are given in [1].
%
%Note that when given only r-u-v and not r-u-v-w coordinates, an implicit
%assumption is that the measurement is not extremely close to being outside
%of the receiver's viewing area. For many practical radars having around a
%+/-60 degree field of view and angular measurement standard deviations
%less than one degress (once u-v deviations are converted to angles), then
%this is not a problem at all. However, given a measurement very close to
%+/- 90 degrees from the boresight, when cubature points are added for the
%integration, this can push u^2+v^2 above 1. At such extremes, however, the
%Gaussian approximation tends to be invalid anyway. Of course, utilizing
%r-u-v-w coordinates (with an often singular covariance matrix) can cause
%this issue to be avoided.
%
%EXAMPLE:
%This example goes over confusion that can arise when using the
%"useHalfRange" option. When useHalfRange=false (the default; get a
%round-trip range), then the noise model is 
% zRT=zTrue+SRrt*randn(3,1);
%where SRrt is SR when considering the round-trip range. That is, the
%measurement is the true round-trip r-u-v value plus the zero-mean Gaussian
%noise. Now, suppose that we were to halve the range. This is the same as
%multipling z by
% C=[1/2,0,0;
%      0,1,0;
%      0,0,1];
%Thus, we get 
% C*z=C*zTrue+C*SRrt*randn(3,1);
%When useHalfRange is true, the assumption is that SR is the covariance
%matrix incorporating C. That is, the measurement model is 
% zMono=zTrue+SRmono*randn(3,1);
%where the one-way square root lower-triangular noise covariance matrix is 
%SRmono=C*SRrt.
%To demonstrate this, we consider a monostatic example, where we use a
%round-trip range and also where we use a one-way range.
% SRrt=diag([10;1e-3;1e-3]);
% zRT=[100e3;0;0];
% %Range is just halved; the measurement is otherwise the same.
% zMono=[100e3/2;0;0];
% [zCartBi,RCartBi]=ruv2CartCubature(zRT, SRrt, false);
% [zCartMonoWrong,RCartMonoWrong]=ruv2CartCubature(zMono, SRrt, true);
% C=[1/2,0,0;
%      0,1,0;
%      0,0,1];
% SRmono=C*SRrt;
% [zCartMono,RCartMono]=ruv2CartCubature(zMono, SRmono, true);
% RCartBi-RCartMonoWrong
% RCartBi-RCartMono
%Thus, one will see that the conversion using the round-trip range only has
%an equivalent Cartesian covariance matrix when one properly adjusts the
%range measurement covariance matrix for the halving.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(xi))
    [xi,w]=fifthOrderCubPoints(3);
end

if(nargin<6||isempty(M))
    M=eye(3);
end

if(nargin<5||isempty(zRx))
    zRx=zeros(3,1);
end

if(nargin<4||isempty(zTx))
    zTx=zeros(3,1);
end

if(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

numMeas=size(z,2);

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,numMeas]);
end

h=@(x)ruv2Cart(x,useHalfRange,zTx,zRx,M);

zCart=zeros(3,numMeas);
RCart=zeros(3,3,numMeas);
for curMeas=1:numMeas
    [zCart(:,curMeas),RCart(:,:,curMeas)]=calcCubPointMoments(z(:,curMeas),SR(:,:,curMeas),h,xi,w);
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
