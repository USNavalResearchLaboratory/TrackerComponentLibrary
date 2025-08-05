function u=polAng2U2D(azimuth,systemType,MP,MU,includeV)
%%POLANG2U2D Convert 2D polar angles into equivalent direction cosine
%            values. The "direction cosine" u is just the x coordinate of a
%            unit vector from the receiver to the target in the coordinate
%            system at the receiver. Optionally, a full [u;v] unit vector
%            in 2D can be returned if includeV is true.
%
%INPUTS: azimuth A 1XN or NX1 set of azimuthal values.
%     systemType An optional parameter specifying the axis from which the
%                angles are measured. Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  The azimuth angle is counterclockwise from the x axis.
%                1 The azimuth angle is measured clockwise from the y axis.
%             MP Optionally, the rotation matrix that would rotate a vector
%                in global coordinates into the local coordinate system
%                defining polar coordinates. If omitted or an empty matrix
%                is passed, it is assumed there is no rotation.
%             MU Optionally, the rotation matrix that would rotate a vector
%                in global coordinates into the local coordinate system
%                into which the u coordinate is computed. If omitted or an
%                empty matrix is passed, it is assumed there is no
%                rotation.
%       includeV An optional boolean value indicating whether a second
%                direction cosine component should be included. The u
%                direction cosine is one parts of a 2D unit vector.
%                Generally, one might assume that the target is in front of
%                the sensor, so the second component would be positive and
%                is not needed. However, the second component can be
%                included if ambiguity exists. The default if this
%                parameter is omitted or an empty matrix is passed is
%                false.
%
%OUTPUTS: u A 1XN (without v) or 2XN (with v) set of direction cosines in
%           2D corresponding to the specified azimuthal values.
%
%EXAMPLE 1:
%This shows that the value of the converted component is consistent with
%what one would get when considering a Cartesian location in both
%coordinate systems.
% %This angle specifies the rotation of the polar coordinate system
% thetaP=2*pi*rand(1);
% MP=rotMat2D(thetaP);
% %This angle specifies the rotation of the u coordinate system.
% thetaU=2*pi*rand(1);
% MU=rotMat2D(thetaU);
% systemType=0;
% thetaPol=deg2rad(35);%The polar angle
% zRx=[100;40];
% xCart=pol2Cart([1;thetaPol],systemType,true,zRx,zRx,MP);
% includeV=true;
% u=getUDirection2D(xCart,zRx,MU,includeV);
% uAlt=polAng2U2D(thetaPol,systemType,MP,MU,includeV);
% RelErr=(uAlt-u)./u
%
%EXAMPLE 2:
%This is the same type of consistency check as the first example, but just
%done in a different manner. The relative error is on the order of finite
%precision limitations.
% tarLoc=[13e3;26e3];
% lRx=[12e3;4e3];
% systemType=0;
% includeV=true;
% Mu=randRotMat(2);
% MAz=randRotMat(2);
% uv=getUDirection2D(tarLoc,lRx,Mu,includeV);
% AzVal=getPolAngle(tarLoc,systemType,lRx,MAz);
% uvConv=polAng2U2D(AzVal,systemType,MAz,Mu,includeV);
% RelErr=max(abs((uvConv-uv)./uv))
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(includeV))
   includeV=false; 
end

if(nargin<4||isempty(MU))
    thetaU=0;
else
    thetaU=rotMat2D2Angle(MU);
end

if(nargin<3||isempty(MP))
    thetaP=0;
else
    thetaP=rotMat2D2Angle(MP);
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Make a row vector.
azimuth=azimuth(:).';

switch(systemType)
    case 0
        azimuth=azimuth-thetaP+thetaU;

        if(includeV)
            u=[cos(azimuth);sin(azimuth)];
        else
            u=cos(azimuth);
        end
    case 1
        azimuth=azimuth+thetaP-thetaU;

        if(includeV)
            u=[sin(azimuth);cos(azimuth)];
        else
            u=sin(azimuth);
        end
    otherwise
        error('Invalid system type specified.')
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
