function azEl=getSpherAngle(zC,systemType,zRx,M)
%%GETSPHERANGLE Convert a position into local spherical angles
%         [azimuth;elevation] in radians ignoring atmospheric and other
%         propagation effects.
%
%INPUT: zC A 3XN matrix of Cartesian points in global [x;y;z] Cartesian
%          coordinates.
% systemType An optional parameter specifying the axes from which the
%          angles are measured in radians. Possible values are
%          0 (The default if omitted) Azimuth is measured
%            counterclockwise from the x-axis in the x-y plane. Elevation
%            is measured up from the x-y plane (towards the z-axis). This
%            is consistent with common spherical coordinate systems for
%            specifying longitude (azimuth) and geocentric latitude
%            (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given
%            elevation, one desires the angle away from the z-axis, which
%            is (pi/2-elevation).
%          3 This is the same as 0 except azimuth is measured clockwise
%            from the y-axis in the x-y plane instead of counterclockwise
%            from the x-axis. This coordinate system often arises when
%            given "bearings" in a local East-North-Up coordinate system,
%            where the bearing directions are measured East of North.
%      zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix is
%          passed, then the receivers are assumed to be at the origin. If
%          only a single vector is passed, then the receiver location is
%          assumed the same for all of the target states being converted.
%        M A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The z-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(3) --the
%          identity matrix is used. If only a single 3X3 matrix is passed,
%          then is=t is assumed to be the same for all N conversions.
%
%OUTPUT: azEl The 2XN matrix of azimuth and elevation angles in radians in
%             the form [azimuth;elevation].
%
%The conversion from Cartesian to spherical coordinates is given in Chapter
%14.4.4.2 of [1].
%
%REFERENCES:
%[1] R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%    to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%    Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(zC,2);

if(nargin<4||isempty(M))
    M=repmat(eye(3),[1,1,N]);
elseif(size(M,3)==1)
    M=repmat(M,[1,1,N]);
end

if(nargin<3||isempty(zRx))
    zRx=zeros(3,N);
elseif(size(zRx,2)==1)
    zRx=repmat(zRx,[1,N]);
end

if(nargin<2||isempty(systemType))
	systemType=0; 
end

azEl=zeros(2,N);
for curPoint=1:N
%The target location in the receiver's coordinate system.
    zCL=M(:,:,curPoint)*(zC(:,curPoint)-zRx(:,curPoint));

%Perform the conversion.
    r1=norm(zCL);%Receiver to target.
    
    x=zCL(1);
    y=zCL(2);
    z=zCL(3);
    
    switch(systemType)
        case 0
            azimuth=atan2(y,x);
            elevation=asin(z./r1);
        case 1
            azimuth=atan2(x,z);
            elevation=asin(y./r1);
        case 2
            azimuth=atan2(y,x);
            elevation=pi/2-asin(z./r1);
        case 3
            azimuth=atan2(x,y);
            elevation=asin(z./r1);
        otherwise
            error('Invalid system type specified')
    end

    azEl(:,curPoint)=[azimuth;elevation];
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
