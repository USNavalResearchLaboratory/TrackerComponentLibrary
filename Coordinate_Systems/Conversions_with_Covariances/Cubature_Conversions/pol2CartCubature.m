function [zCart,RCart]=pol2CartCubature(z,SR,systemType,useHalfRange,zTx,zRx,M,xi,w)
%%POL2CARTCUBATURE Use cubature integration to approximate the moments of a
%          noise-corrupted bistatic polar location converted into 2D
%          Cartesian coordinates. This function ignores all propagation
%          effects, such as atmospheric refraction. The angle can measured
%          either counterclockwise from the x-axis, which is standard in
%          mathematics, or clockwise from the y axis, which is more common
%          in navigation.
%
%INPUTS: z A 2XN set of N polar points of the form [range;azimuth] that are
%          to be converted into Cartesian coordinates. The azimuth is in
%          radians.
%       SR The 2X2XN lower-triangular square roots of the covariance
%          matrices associated with z. If all of the matrices are the same,
%          then this can just be a single 2X2 matrix.
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are
%          0 (The default if omitted or an empty matrix is passed) The
%            azimuth angle is counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided, or an
%          empty matrix is passed, is true.
%      zTx The 2X1 [x;y] location vector of the transmitter in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitter is assumed to be at the
%          origin. zTx can have more than 3 rows; additional rows are
%          ignored.
%      zRx The 2X1 [x;y] location vector of the receiver in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix is
%          passed, then the receiver is assumed to be at the origin. zRx
%          can have more than 3 rows; additional rows are ignored.
%        M A 2X2 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted or an
%          empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used.
%       xi A 2XnumCubaturePoints (for position-only) matrix of cubature
%          points for the numeric integration. If this and the final
%          parameter are omitted or empty matrices are passed, then
%          fifthOrderCubPoints is used to generate cubature points.
%        w A numCubaturePointsX1 vector of the weights associated with the
%          cubature points.
%
%OUTPUTS: zCart The 2XN approximate means of the PDF of the polar
%               measurements converted to [x;y] Cartesian coordinates.
%         RCart The 2X2XN set of approximate covariance matrices for the N
%               estimates.
%
%Details of the numerical integration used in the conversion are given in
%[1].
%
%EXAMPLE:
%This example goes over confusion that can arise when using the
%"useHalfRange" option. When useHalfRange=false, then the noise model is 
% zRT=zTrue+SRrt*randn(2,1);
%where SRrt is SR when considering the round-trip range. That is, the
%measurement is the true round-trip r-theta value plus the zero-mean
%Gaussian noise. Now, suppose that we were to halve the range. This is the
%same as multipling z by
% C=[1/2,0;
%      0,1];
%Thus, we get 
% C*z=C*zTrue+C*SRrt*randn(2,1);
%When useHalfRange is true, the assumption is that SR is the covariance
%matrix incorporating C. That is, the measurement model is 
% zMono=zTrue+SRmono*randn(3,1);
%where the one-way square root lower-triangular noise covariance matrix is 
%SRmono=C*SRrt.
%To demonstrate this, we consider a monostatic example, where we use a
%round-trip range and also where we use a one-way range.
% SRrt=diag([10;1e-3]);
% zRT=[100e3;0];
% %Range is just halved; the measurement is otherwise the same.
% zMono=[100e3/2;0];
% [zCartBi,RCartBi]=pol2CartCubature(zRT, SRrt,0,false);
% [zCartMonoWrong,RCartMonoWrong]=pol2CartCubature(zMono, SRrt, true);
% C=[1/2,0;
%      0,1];
% SRmono=C*SRrt;
% [zCartMono,RCartMono]=pol2CartCubature(zMono, SRmono,0,true);
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
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(z,2);
    
    if(nargin<8||isempty(xi))
        [xi,w]=fifthOrderCubPoints(2);
    end

    if(nargin<7||isempty(M))
        M=eye(2);
    end

    if(nargin<6||isempty(zRx))
        zRx=zeros(2,1);
    end

    if(nargin<5||isempty(zTx))
        zTx=zeros(2,1);
    end

    if(nargin<4||isempty(useHalfRange))
        useHalfRange=true;
    end

    if(nargin<3||isempty(systemType))
        systemType=0; 
    end
    
    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numPoints]);
    end
    
    h=@(z)pol2Cart(z,systemType,useHalfRange,zTx(1:2,:),zRx(1:2,:),M);
    
    zCart=zeros(2,numPoints);
    RCart=zeros(2,2,numPoints);
    for curMeas=1:numPoints
        [zCart(:,curMeas), RCart(:,:,curMeas)]=calcCubPointMoments(z(:,curMeas),SR(:,:,curMeas),h,xi,w);
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
