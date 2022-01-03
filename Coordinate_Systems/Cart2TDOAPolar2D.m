function z=Cart2TDOAPolar2D(points,systemType,lRef,lRx,M,c)
%%CART2TDOAPOLAR2D Given 2DCartesion positions, convert them to a
%       time-difference of arrival (TDOA) component and a polar angle
%       The angle is assumed measured at the non-reference receiver.
%
%INPUTS: points A 2XN set of Cartesian points.
%  systemType An optional parameter specifying the axes from which the
%           angles are measured. Possible values are
%           0 (The default if omitted) The azimuth angle is
%             counterclockwise from the x axis.
%           1 The azimuth angle is measured clockwise from the y axis.
%      lRef The 2X1 location of the reference receiver in global
%           coordinates.
%       lRx The 2X1 location of the receiver in global coordinates at which
%           the direction of arrival in spherical coordinates is locally
%           defined and whose time of arrival (TOA) has the TOA at lRef
%           subtracted to get the TDOA.
%         M A 2X2 rotation matrix to go from the alignment of the global
%           coordinate system to that at the receiver. If omitted, then it
%           is assumed that the local coordinate system is aligned with the
%           global and M=eye(2) --the identity matrix is used.
%         c The propagation speed in the medium in question. If this
%           parameter is omitted or an empty matrix is passed, the default
%           value of Constants.speedOfLight is used.
%
%OUTPUTS: z A 2XN set of points of the form [TDOA;azimuth;elevation] with
%           the angles given in radians.
%
%This function just calls getTDOA and getPolAngle
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<5||isempty(M))
    M=eye(2,2);
end

if(isempty(systemType))
    systemType=0;
end

z=[getTDOA(points,lRef,lRx,c);
   getPolAngle(points,systemType,lRx,M)];

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
