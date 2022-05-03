function [zCart,RCart]=spher2CartCubature(spherPoint,SR,systemType,useHalfRange,zTx,zRx,M,xi,w)
%%SPHER2CARTCUBATURE Use cubature integration to approximate the moments
%                    of a noise-corrupted spherical location converted into
%                    3D Cartesian coordinates. An option allows for
%                    the angles to be specified from different axes.
%                    Optionally, a bistatic range can be used when
%                    considering bistatic measurements in a local spherical
%                    coordinate system.
%
%INPUTS: spherPoint One or more points given in terms of range, azimuth and
%           elevation, with the angles in radians To convert N points,
%           spherPoint is a 3XN matrix with each column having the
%           format [range;azimuth; elevation].
%        SR The 3X3XN lower-triangular square roots of the covariance
%           matrices associated with polPoint. If all of the matrices
%           are the same, then this can just be a single 3X3 matrix.
% systemType An optional parameter specifying the axis from which the
%           angles are measured. Possible values are
%           0 (The default if omitted) Azimuth is measured
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
% useHalfRange An optional boolean value specifying whether the bistatic
%           (round-trip) range value has been divided by two. This normally
%           comes up when operating in monostatic mode (the most common
%           type of spherical coordinate system), so that the range
%           reported is a one-way range (or just half a bistatic range).
%           The default if this parameter is not provided is false if zTx
%           and zRx are provided and true if they are omitted (monostatic).
%       zTx The 3XN [x;y;z] location vectors of the transmitters in global
%           Cartesian coordinates. If this parameter is omitted, then the
%           transmitters are assumed to be at the origin. If only a single
%           vector is passed, then the transmitter location is assumed the
%           same for all of the target states being converted. zTx can have
%           more than 3 rows; additional rows are ignored. If monostatic
%           an empty matrix can be passed.
%       zRx The 3XN [x;y;z] location vector of the receivers in Cartesian
%           coordinates.  If this parameter is omitted, then the
%           receivers are assumed to be at the origin. If only a single
%           vector is passed, then the receiver location is assumed the
%           same for all of the target states being converted. zRx can have
%           more than 3 rows; additional rows are ignored. If monostatic
%           an empty matrix can be passed.
%        M  A 3X3XN hypermatrix of the rotation matrices to go from the
%           alignment of the global coordinate system to that at the
%           receiver. The z-axis of the local coordinate system of the
%           receiver is the pointing direction of the receiver. If omitted,
%           then it is assumed that the local coordinate system is aligned
%           with the global and M=eye(3) --the identity matrix is used. If
%           only a single 3X3 matrix is passed, then it is assumed to be
%           the same for all of the N conversions.
%        xi A 3 X numCubaturePoints matrix of cubature points for the
%           numeric integration. If this and the next parameter are
%           omitted or empty matrices are passed, then
%           fifthOrderCubPoints is used to generate cubature points.
%         w A numCubaturePoints X 1 vector of the weights associated
%           with the cubature points.
%
%OUTPUTS: zCart The approximate means of the PDF of the Cartesian converted
%               measurements in [x;y;z] Cartesian coordinates for each
%               measurement. This is a 3XN matrix.
%         RCart The approximate 3X3XN set of covariance matrices of the
%               PDFs of the Cartesian converted measurements.
%
%Details of the numerical integration used in the conversion are given in
%[1].
%
%EXAMPLE:
%This example goes over confusion that can arise when using the
%"useHalfRange" option. When useHalfRange=false, then the noise model is 
% zRT=zTrue+SRrt*randn(3,1);
%where SRrt is SR when considering the round-trip range. That is, the
%measurement is the true round-trip r-azimuth-elevation value plus the
%zero-mean Gaussiannoise. Now, suppose that we were to halve the range.
%This is the same as multipling z by
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
% [zCartBi,RCartBi]=spher2CartCubature(zRT, SRrt,0,false);
% [zCartMonoWrong,RCartMonoWrong]=spher2CartCubature(zMono, SRrt,0,true);
% C=[1/2,0,0;
%      0,1,0;
%      0,0,1];
% SRmono=C*SRrt;
% [zCartMono,RCartMono]=spher2CartCubature(zMono, SRmono,0,true);
% RCartBi-RCartMonoWrong
% RCartBi-RCartMono
%Thus, one will see that the conversion using the round-trip range only has
%an equivalent Cartesian covariance matrix when one properly adjusts the
%range measurement covariance matrix for the halving.
%
%REFERENCES:
%[1] David F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(systemType))
        systemType=0;
    end

    numDim=size(spherPoint,1);
    numPoints=size(spherPoint,2);

    if(numDim~=3)
        error('The spherical points have the wrong dimensionality')
    end

    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numPoints]);
    end
    
    if(nargin<8||isempty(xi))
       [xi,w]=fifthOrderCubPoints(numDim);
    end
    
    if(nargin<5||isempty(zTx))
        zTx=[];
    end

    if(nargin<6||isempty(zRx))
        zRx=[];
    end

    if(nargin<7||isempty(M))
        M=[]; 
    end
    
    if(isempty(zTx)&&(nargin<4||isempty(useHalfRange)))
        useHalfRange=true;
    elseif(~isempty(zTx)&&(nargin<4||isempty(useHalfRange)))
        useHalfRange=false;
    end
    
    h=@(x)spher2Cart(x,systemType,useHalfRange,zTx,zRx,M);
    
    %Allocate space for the return variables.
    zCart=zeros(numDim,numPoints);
    RCart=zeros(numDim,numDim,numPoints);
    for curPoint=1:numPoints
        [zCart(:,curPoint),RCart(:,:,curPoint)]=calcCubPointMoments(spherPoint(:,curPoint),SR(:,:,curPoint),h,xi,w);
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
