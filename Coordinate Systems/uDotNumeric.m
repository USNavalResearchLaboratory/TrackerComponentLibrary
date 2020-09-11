function uDot=uDotNumeric(u,x,t,vertFunc,epsVal)
%%UDOTNUMERIC Use numerical differentiation to obtain the derivative of the
%             basis vectors for the local cordinate system of a moving
%             observer on a curved surface. Such an evolving coordinate
%             system has been termed "wander coordinate" or "naturally
%             evolving coordinates". Targets moving in non-maneuvering,
%             level flight with a coordinate system evolving as per this
%             coordinate system will follow geodesic curves. When traveling
%             on/ over a reference ellipsoid, the function uDotEllipsoid is
%             preferred as the solution is available without numeric
%             integration.
%
%INPUTS: u The orthonormal basis vectors for the local coordinate system.
%          u(:,3) is the local vertical.
%        x An nX1 vector whose first three elements are Cartesian position
%          in the global coordinate system and whose next three elements
%          are velocity in the local coordinate system. Other elements of x
%          do not matter.
%        t A time component that is passed to vertFunc, in case the local
%          vertical is somehow time-dependent.
% downFunc A function handle of the format vertFunc(r,t), where r is a
%          location in global coordinates and t is the above time
%          parameter. The function handle returns a vector pointing in the
%          direction of the local "down" direction. The vector does not
%          need to be unit magnitude. For example, an acceleration due to
%          gravity could be passed.
%   epsVal An optional parameter specifying the step size taken in each
%          dimension for numerical differentiation. If this parameter is
%          omitted, then a step size of max(1e-4*norm(x(1:3)),1e-10) is
%          used, which is reasonable if x has a position on the reference
%          ellipsoid.
%
%OUTPUTS: uDot The derivative of the basis vector u with respect to time.
%
%The solution is detailed in [1]. Compare to uDotEllipsoid where an
%explicit solution for an ellipsoidal surface is possible.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(epsVal))
        epsVal=max(1e-4*norm(x(1:3)),1e-10);
    end

    r=x(1:3);%The location in global coordinates.

    %The local vertical
    g=vertFunc(r);    

%The Jacobian of the local vertical
    theVert=@(r)vertFunc(r,t);
    H=numDiff(r,theVert,3,3,epsVal);
    
%Perform the appropriate transformations of the Jacobian of the local
%vertical to obtain the time-derivatives of the local coordinate axes.

    %Convert the local velocity into a global velocity.
    velGlob=getGlobalVectors(x(4:6),u);
    gDot=H*velGlob;
    normG=norm(g);
    
    A=[-u(:,2),u(:,1)];
    
    c=lsqminnorm(A,(-gDot/normG+(g/normG^3)*dot(gDot,g)));
    Omega=c(1)*u(:,1)+c(2)*u(:,2);
    
    uDot(:,1)=cross(Omega,u(:,1));
    uDot(:,2)=cross(Omega,u(:,2));
    uDot(:,3)=cross(Omega,u(:,3));
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
