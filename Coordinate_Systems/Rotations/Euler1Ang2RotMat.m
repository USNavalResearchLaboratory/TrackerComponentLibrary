function M=Euler1Ang2RotMat(theta,series,handed)
%%EULER1ANG2ROTMAT  Obtain the Euler-angle rotation matrix for a rotation
%            of angle theta about single axis (x, y, or z). The rotation is
%            either right-handed (counterclockwise, the default) or
%            left-handed depending on handed. The components of vectors to
%            be rotated are assumed ordered [x;y;z].
%
%INPUTS: theta The angle in radians about which the rotation matrix should
%              rotate a vector around the axis given by series. The
%              handedness of the angle is given by handed.
%       series A text string specifying the angle about which the rotation
%              is performed. This can be 'x', 'y', or 'z'.
%       handed The handedness of the rotation angle. If omitted, it is
%              assumed that the rotation is right-handed (the standard).
%              Possible values are:
%              'right' The default if omitted. The rotation is right-
%                      handed.
%              'left'  The rotation is left-handed. The rotation angle is
%                      clockwise when one is looking into the rotation
%                      axis.
%
%OUTPUTS: M The rotation matrix such that M*v rotates a vector v by theta
%           about the specified axis.
%
%Euler angles are discussed in [1].
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. - Dec. 1993.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(handed))
   handed='right';
end

switch(handed)
    case 'right'
    case 'left'
        theta=-theta;
    otherwise
        error('Invalid handedness provided.')
end

cosT=cos(theta);
sinT=sin(theta);

switch(series)
    case 'x'
        M=[1, 0,       0;
           0, cosT, -sinT;
           0, sinT,  cosT];
    case 'y'
        M=[cosT, 0,  sinT;
              0, 1,     0;
          -sinT, 0, cosT];
    case 'z'
        M=[cosT, -sinT, 0;
           sinT,  cosT, 0;
              0,    0, 1];
    otherwise
        error('Invalid rotation axis specified')
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
