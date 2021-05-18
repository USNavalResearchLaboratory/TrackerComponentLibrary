function M=Euler2Ang2RotMat(theta1,theta2,series,handed)
%%EULER2ANG2ROTMAT Get a rotation matrix for a series of 2 Euler angles.
%              The axes of rotation are given by series. The rotations are
%              in order of theta2, and then theta1 about the specified
%              axes. Note that the axis of rotation after theta2 is a
%              ROTATED axes. The rotations are either right or left-handed,
%              with right-handed rotations being the default. The
%              components of vectors to be rotated are assumed ordered
%              [x;y;z].
%
%INPUTS: theta,theta2 The two angles of rotation given in radians about the
%              two axes specified in series the rotations are taken. The
%              first rotation is theta2, then theta1.
%       series A character string specifying the series of axes about which
%              the two rotations are taken. All possible combinations of
%              axes without repeating an axis are valid. For example, 'xy'
%              means rotate theta2 about the y axis, then rotate theta1
%              about the rotated x axis. All possible combinations of values
%              are: 'xy', 'xz', 'yx', 'yz', 'zx', and 'zy'.
%       handed The handedness of the rotation angle. If omitted, it is
%              assumed that the rotation is right-handed (the standard).
%              Possible values are:
%              'right' The default if omitted. The rotation is right-
%                      handed.
%              'left'  The rotation is left-handed. The rotation angle is
%                      clockwise when one is looking into the rotation
%                      axis.
%
%OUTPUTS: M The 3X3 rotation matrix such that M*v rotates the 3X1 vector v
%           according to the given rotation angles and axes.
%
%Euler angles are discussed in [1].
%
%EXAMPLE:
%Here, we rotate by some arbitrary angles and one can see that the
%rotations from this functions are equal to sunsequently appending the
%rotations.
% theta1=2*pi*rand();
% theta2=2*pi*rand();
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'xy','right')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'y'))))
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'xz','right')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'z'))))
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'yx','right')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'x'))))
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'yz','right')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'z'))))
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'zx','right')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'x'))))
% max(max(abs(Euler2Ang2RotMat(theta1,theta2,'zy','right')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'y'))))
%All of the above differences will be zero or around eps(), indicating
%agreement within finite precision limits.
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. - Dec. 1993.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(handed))
    handed='right';
end

switch(handed)
    case 'right'
    case 'left'
        theta1=-theta1;
        theta2=-theta2;
    otherwise
        error('Invalid handedness provided.')
end

c1=cos(theta1);
s1=sin(theta1);
c2=cos(theta2);
s2=sin(theta2);

switch(series)
    case 'xy'
        M=[    c2,  0,     s2;
            s1*s2, c1, -c2*s1;
           -c1*s2, s1,  c1*c2];
    case 'xz'
        M=[   c2,   -s2,   0;
           c1*s2, c1*c2, -s1;
           s1*s2, c2*s1,  c1];
    case 'yx'
        M=[ c1, s1*s2, c2*s1;
             0,    c2,   -s2;
           -s1, c1*s2, c1*c2];
    case 'yz'
        M=[ c1*c2, -c1*s2, s1;
               s2,     c2,  0;
           -c2*s1,  s1*s2, c1];
    case 'zx'
        M=[c1, -c2*s1,  s1*s2;
           s1,  c1*c2, -c1*s2;
            0,     s2,     c2];
    case 'zy'
        M=[c1*c2, -s1, c1*s2;
           c2*s1,  c1, s1*s2;
             -s2,   0,    c2];
    otherwise
        error('Invalid rotation series provided.')
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
