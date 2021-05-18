function vRot=rotateVectorAxisAng(v,u,theta,handed)
%%ROTATEVECTORAXISANG Given three-dimensionals vector v and a unit vector
%                     u, rotate the vectors an angle of theta
%                     counterclockwise (right-handed) or clockwise
%                     (left-handed) about the vector (axis) u.
%
%INPUTS: v The 3XN set of N vectors that are to be rotated.
%        u A unit vector representing the axis about which v is to be
%          rotated.
%    theta The angle in radians by which v is to be rotated about u. When
%          right-handed, the rotation angle is clockwise when one is
%          looking in the same direction that the rotation axis points. A
%          left-handed rotation is the opposite.
%   handed The handedness of the rotation angle. If omitted, it is assumed
%          that the rotation is right-handed (the standard). Possible
%          values are
%  'right' The default if omitted. The rotation is right-handed.
%   'left' The rotation is left-handed.
%
%OUTPUTS: vRot The 3XN set of vectors v rotated an angle of theta about the
%              axis u according to the given handedness.
%
%This simply implements the Rodrigues' rotation formula. The formula is
%given in equations (96) and (97) in [1], where changes to allow for both
%left and right-handed rotations have been made.
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%the Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct.-Dec. 1993.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(handed))
   handed='right'; 
end

switch(handed)
    case 'right'
    case 'left'
        theta=-theta;
    otherwise
        error('Invalid handedness provided.')
end

N=size(v,2);
vRot=zeros(3,N);%Allocate space

for curVec=1:N
    vRot(:,curVec)=v(:,curVec)*cos(theta)+cross(u,v(:,curVec))*sin(theta)+u*dot(u,v(:,curVec))*(1-cos(theta));
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
