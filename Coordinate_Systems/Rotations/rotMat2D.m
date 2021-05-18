function M=rotMat2D(theta,isClockwise)
%%ROTMAT2D Obtain a rotation matrix in two dimensions (x-y Cartesian space)
%          for a given angle.
%
%INPUTS: theta The angle in radians about which the rotation matrix should
%              rotate a 2X1 vector. 
%  isClockwise An optional boolean parameter indication whether the
%              rotation is clockwise. The default if omitted or an empty
%              matrix is passed is false.
%
%OUTPUTS: M The 2X2 rotation matrix.
%
%The expression for the 2D rotation matrix is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Rotation Matrix." From MathWorld--A Wolfram Web
%   Resource. http://mathworld.wolfram.com/RotationMatrix.html
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(isClockwise))
    isClockwise=false; 
end

if(isClockwise)
    theta=-theta;
end

cosTheta=cos(theta);
sinTheta=sin(theta);

M=[cosTheta, -sinTheta;
   sinTheta, cosTheta];
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
