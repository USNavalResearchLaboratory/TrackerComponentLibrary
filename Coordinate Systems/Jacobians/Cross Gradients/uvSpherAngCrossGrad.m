function JTotal=uvSpherAngCrossGrad(azEl,systemType,includeW)
%%UVSPHERANGCROSSGRAD Determine the partial derivatives of 3D spherical
%                  angular components with respect to u-v direction
%                  cosines.
%
%INPUTS: azEl A 2XN set of points in 3D [azimuth;elevation] in radians.
% systemType An optional parameter specifying the axes from which the
%          angles for the spherical coordinate system are measured in
%          radians. Possible vaues are
%          0 (The default if omitted) Azimuth is measured counterclockwise
%            from the x-axis in the x-y plane. Elevation is measured up
%            from the x-y plane (towards the z-axis). This is consistent
%            with common spherical coordinate systems for specifying
%            longitude (azimuth) and geocentric latitude (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given
%            elevation, one desires the angle away from the z-axis, which
%            is (pi/2-elevation).
% includeW If one wants to have partial derivatives with respect to the
%          third component of the direction cosine unit vector (the z-axis
%          component), then set this to true. the default if omitted or an
%          empty matrix is passed is false.
%
%OUTPUTS: JTotal A 2X2XnumPoints (or 3X2XnumPoints) set of Jacobian
%            matrices where the rows are u and v (and w if included) and
%            the columns take the derivative of the row components with
%            respect to azimuth and v in that order.
%
%The Jacobian is a matrix of partial derivatives. The Jacobian returned by
%spherAngUvCrossGrad is the matrix inverse of that returned by this
%function.
%
%EXAMPLE:
%Here, we verify that the partial derivatives returned by this function are
%about equal to those returned via numeric differentiation (forward
%differencing).
% points=[[0.1;0.2],[0.1;0.1],[0;0]];%AzEl points
% epsVal=1e-9;
% systemType=2;
% 
% uv=spherAng2Uv(points,systemType);
% uv1=spherAng2Uv(bsxfun(@plus,points,[epsVal;0]),systemType);
% numDiffAz=(uv1-uv)/epsVal;
% uv1=spherAng2Uv(bsxfun(@plus,points,[0;epsVal]),systemType);
% numDiffEl=(uv1-uv)/epsVal;
% gradVal=uvSpherAngCrossGrad(points,systemType);
% 
% max(abs(numDiffAz(:)-vec(gradVal(:,1,:))))
% max(abs(numDiffEl(:)-vec(gradVal(:,2,:))))
%One will see that both numeric differences are on the order of 9e-9 and
%1.5e-7, which is a good amount of agreement.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
   systemType=0; 
end

if(nargin<3||isempty(includeW))
    includeW=false;
end

N=size(azEl,2);

if(includeW)
    JDim=3;
else
    JDim=2;
end

JTotal=zeros(JDim,2,N);
for curPoint=1:N
    az=azEl(1,curPoint);
    el=azEl(2,curPoint);
    
    sinAz=sin(az);
    cosAz=cos(az);
    sinEl=sin(el);
    cosEl=cos(el);
    
    J=zeros(JDim,2);
    switch(systemType)
        case 0
            %Azimuth derivatives
            J(1,1)=-cosEl*sinAz;
            J(2,1)=cosAz*cosEl;

            %Elevation derivatives
            J(1,2)=-cosAz*sinEl;
            J(2,2)=-sinAz*sinEl;

            if(includeW)
                J(3,1)=0;
                J(3,2)=cosEl;
            end
        case 1
            %Azimuth derivatives
            J(1,1)=cosAz*cosEl;
            J(2,1)=0;

            %Elevation derivatives
            J(1,2)=-sinAz*sinEl;
            J(2,2)=cosEl;

            if(includeW)
                J(3,1)=-cosEl*sinAz;
                J(3,2)=-cosAz*sinEl;
            end
        case 2
            %Azimuth derivatives
            J(1,1)=-sinAz*sinEl;
            J(2,1)=cosAz*sinEl;

            %Elevation derivatives
            J(1,2)=cosAz*cosEl;
            J(2,2)=cosEl*sinAz;

            if(includeW)
                J(3,1)=0;
                J(3,2)=-sinEl;
            end
        otherwise
            error('Invalid system type specified.')
    end
    JTotal(:,:,curPoint)=J;
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
