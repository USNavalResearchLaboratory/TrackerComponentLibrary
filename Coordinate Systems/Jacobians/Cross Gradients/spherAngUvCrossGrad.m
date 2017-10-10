function JTotal=spherAngUvCrossGrad(uv,systemType)
%%SPHERANGUVCROSSGRAD Determine the partial derivatives of u-v direction
%                  cosine components with respect to 3D spherical angles.
%
%INPUTS: uv A 2XnumPoints or 3XnumPoints (if the third components of the
%          unit vector) set of direction [u;v;w] cosines values in 3D. If
%          the third component of the unit vector is omitted, it is assumed
%          to be positive.
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
%
%OUTPUTS: JTotal A 2X2XnumPoints set of Jacobian matrices where the rows
%            are azimuth and elevation the columns take the derivative of
%            the row components with respect to u and v in that order. Note
%            that singularities exist, in which case NaNs can be returned.
%
%The Jacobian is a matrix of partial derivatives. The Jacobian returned by
%uvSpherAngCrossGrad is the matrix inverse of that returned by this
%function.
%
%EXAMPLE:
%Here, we verify that the partial derivatives returned by this function are
%about equal to those returned via numeric differentiation (forward
%differencing).
% points=[[0.1;0.2],[0.1;0.1],[0;0]];%AzEl points
% epsVal=1e-9;
% systemType=0;
% 
% azEl=uv2SpherAng(points,systemType);
% azEl1=uv2SpherAng(bsxfun(@plus,points,[epsVal;0]),systemType);
% numDiffU=(azEl1-azEl)/epsVal;
% azEl1=uv2SpherAng(bsxfun(@plus,points,[0;epsVal]),systemType);
% numDiffV=(azEl1-azEl)/epsVal;
% gradVal=spherAngUvCrossGrad(points,systemType);
% 
% max(abs(numDiffU(:)-vec(gradVal(:,1,:))))
% max(abs(numDiffV(:)-vec(gradVal(:,2,:))))
%One will see that both numeric differences are on the order of 2e-7, which
%is a good amount of agreement.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
   systemType=0; 
end

hasW=size(uv,1)>2;

N=size(uv,2);

JTotal=zeros(2,2,N);
for curPoint=1:N
    u=uv(1,curPoint);
    v=uv(2,curPoint);
    
    if(hasW)
        w=uv(3,curPoint);
    else
        w=sqrt(1-u^2-v^2);
    end
    J=zeros(2,2);
    switch(systemType)
        case 0
            %Derivatives with respect to u.
            J(1,1)=-v/(u^2+v^2);
            J(2,1)=-u/(w*sqrt(u^2+v^2));
            %Derivatives with respect to v.
            J(1,2)=u/(u^2+v^2);
            J(2,2)=-v/(w*sqrt(u^2+v^2));
        case 1
            %Derivatives with respect to u.
            J(1,1)=1/w;
            J(2,1)=0;
            %Derivatives with respect to v.
            J(1,2)=u*v/((1-v^2)*w);
            J(2,2)=1/sqrt(1-v^2);
        case 2
            %Derivatives with respect to u.
            J(1,1)=-v/(u^2+v^2);
            J(2,1)=u/(w*sqrt(u^2+v^2));
            %Derivatives with respect to v.
            J(1,2)=u/(u^2+v^2);
            J(2,2)=v/(w*sqrt(u^2+v^2));
        otherwise
            error('Invalid system type specified')
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
