function JTotal=polAngGradient(xG,systemType,lRx,M)
%%POLANGGRADIENT Determine the gradient of the angle of a 2D polar
%           measurement with respect to position (gradient components for
%           velocity etc. are zero and are not provided).Atmospheric and
%           other propagation effects are not taken into account. The
%           measurement can be monostatic or bistatic, but the location of
%           the transmitter does not matter.
%
%INPUTS: xG A 2XN set of Cartesian locations in the format [x;y] where the
%          gradient is desired.
% systemType An optional parameter specifying the axis from which the
%          azimuth angle is measured. It is assumed that the azimuth
%          angle is given in radians. Possible values are
%          0 (The default if omitted) The azimuth angle is
%             counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
%      lRx The 2X1 [x;y] location vector of the receiver in Cartesian
%          coordinates. If this parameter is omitted or an empty matrix is
%          passed, then the receiver is assumed to be at the origin.
%        M A 2X2 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted or an
%          empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2,2)
%          --the identity matrix is used. 
%
%OUTPUTS: JTotal A 1X2XN set of gradient matrices of the polar angle taken
%           with respect to the components [x,y] in 2D in that order for
%           each of the N points in xG.
%
%The conversion utilizing bistatic polar measurements in 2D is similar to
%that using bistatic r-u-v measurements in 3D, which is discussed in [1].
%Basic differentiation was analytically performed to obtain the expressions
%used in this file.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(xG,2);

if(nargin<4||isempty(M))
   M=eye(2,2); 
end

if(nargin<3||isempty(lRx))
   lRx=zeros(2,1);
end

if(nargin<2||isempty(systemType))
   systemType=0; 
end

JTotal=zeros(1,2,N);
for curPoint=1:N
    %Convert the point into the local coordinate system of the receiver.
    xLocal=M*(xG(1:2,curPoint)-lRx(1:2));

    xL=xLocal(1);
    yL=xLocal(2);

    J=zeros(1,2);
    switch(systemType)
        case 0
            %Derivative with respect to x.
            J(1,1)=-yL/(xL^2+yL^2);

            %Derivative with respect to y.
            J(1,2)=xL/(xL^2+yL^2);
        case 1
            %Derivative with respect to x.
            J(1,1)=yL/(xL^2+yL^2);

            %Derivative with respect to y.
            J(1,2)=-xL/(xL^2+yL^2);
        otherwise
            error('Invalid system type specified.')
    end

    %Undo the rotation to move the gradient into the global coordinate
    %system.
    JTotal(:,:,curPoint)=J*M;
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
