function J=uPolar2DCrossGrad(azimuth,systemType,MP,MU,includeV)
%%UPOLAR2DCROSSGRAD Given a 2D azimuthal angle, obtain the derivative of
%              the direction cosine value u with respect to the azimuthal
%              angle. This is the derivative of the output of polAng2U2D
%              with respect to its input.
%
%INPUTS: azimuth A 1XN or NX1 set of N polar angles in radians.
%   systemType An optional parameter specifying the axis from which the
%              angles are measured. Possible values are:
%              0 (The default if omitted or an empty matrix is passed) The
%                azimuth angle is counterclockwise from the x axis.
%              1 The azimuth angle is measured clockwise from the y axis.
%           MP Optionally, the rotation matrix that would rotate a vector
%              in global coordinates into the local coordinate system
%              defining polar coordinates. If omitted or an empty matrix is
%              passed, it is assumed there is no rotation.
%           MU Optionally, the rotation matrix that would rotate a vector
%              in global coordinates into the local coordinate system
%              into which the u cooridnate is computed. If omitted or an
%              empty matrix is passed, it is assumed there is no rotation. 
%     includeV An optional boolean value indicating whether the
%              derivative with respect to a second direction cosine
%              component should be included. The u direction cosine is
%              one parts of a 2D unit vector. If this is true, then the
%              other component is included. The default if omitted or an
%              empty matrix is passed is false.
%
%OUTPUTS: J A 1XnumPoints (if includeV is false) or a 2XnumPoints matrix
%          (if includeV is true) where the rows are u and v differentiated
%          with respect to the azimuthal angle at each point.
%
%EXAMPLE:
%Here, we verify that the partial derivatives returned by this function are
%about equal to those returned via numeric differentiation (forward
%differencing).
% %This angle specifies the rotation of the polar coordinate system
% thetaP=2*pi*rand(1);
% MP=rotMat2D(thetaP);
% %This angle specifies the rotation of the u coordinate system.
% thetaU=2*pi*rand(1);
% MU=rotMat2D(thetaU);
% points=[0.1,-0.2,0,-pi,pi,pi/2];%Az points
% systemType=1;
% epsVal=1e-9;
% includeV=true;
% u=polAng2U2D(points,systemType,MP,MU,includeV);
% u1=polAng2U2D(points+epsVal,systemType,MP,MU,includeV);
% JNumDiff=(u1-u)/epsVal;
% J=uPolar2DCrossGrad(points,systemType,MP,MU,includeV);
% MaxErr=max(max(abs(JNumDiff-J)./J))
%One will see that the disagreement is on the order of 8e-8, which
%indicates a good agreement.
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

if(includeV)
    JDim=2;
else
    JDim=1;
end

N=length(azimuth);

%Make a row vector.
azimuth=azimuth(:).';

J=zeros(JDim,N);
switch(systemType)
    case 0
        azimuth=azimuth-thetaP+thetaU;

        J(1,:)=-sin(azimuth);
        if(includeV)
            J(2,:)=cos(azimuth);
        end
    case 1
        azimuth=azimuth+thetaP-thetaU;

        J(1,:)=cos(azimuth);
        if(includeV)
            J(2,:)=-sin(azimuth);
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
