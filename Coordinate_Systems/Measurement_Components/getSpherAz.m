function az=getSpherAz(zC,systemType,zRx,M)
%%GETSPHERAZ Obtain only the azimuthal component of a position converted
%         into local spherical angles [azimuth;elevation] in radians
%         ignoring atmospheric and other propagation effects.
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
%      zRx The 3X1 [x;y;z] location of the receiver in Cartesian
%          coordinates. This is assumed to be the same for all
%          measurements.  If this parameter is omitted or an empty matrix
%          is passed, then the receiver is assumed to be at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to that at the receiver. The z-axis of the
%          local coordinate system of the receiver is often the pointing
%          direction of the receiver. If omitted or an empty matrix is
%          passed, then it is assumed that the local coordinate system is
%          aligned with the global and M=eye(3).
%
%OUTPUT: az The 1XN vector of azimuth angles in radians.
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
    M=eye(3);
end

if(nargin<3||isempty(zRx))
    zRx=zeros(3,1);
end

if(nargin<2||isempty(systemType))
	systemType=0; 
end

az=zeros(1,N);
for curPoint=1:N
%The target location in the receiver's coordinate system.
    zCL=M*(zC(:,curPoint)-zRx);

%Perform the conversion.
    x=zCL(1);
    y=zCL(2);
    z=zCL(3);
    switch(systemType)
        case 0
            az(curPoint)=atan2(y,x);
        case 1
            az(curPoint)=atan2(x,z);
        case 2
            az(curPoint)=atan2(y,x);
        case 3
            az(curPoint)=atan2(x,y);
        otherwise
            error('Invalid system type specified')
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
