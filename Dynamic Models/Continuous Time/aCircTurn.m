function [aDeriv,aJacob,aHess,papt]=aCircTurn(x,uTurn)
%%ACIRCTURN The drift function for a constant-speed circular turning model
%           about an arbitrary axis.
%
%INPUTS: x The 7X1 target state at time t. It consists of 3D position,
%          velocity, and the scalar turn rate. The axis of the turn is a
%          parameter to this drift function.
%    uTurn A unit vector indicating the axis of rotation about which the
%          turn is performed. If this parameter is omitted or an empty
%          matrix is passed, then uTurn=[0;0;1] is used, which means that
%          the turn takes place in the x-y plane.
%
%OUTPUTS: aDeriv The 7X1 flat-Earth time-derivative of the state.
%       aJacob A 7X7 matrix of first partial derivatives of aVal with respect
%              to the elements of xState. aJacob(:,k) is the derivative of
%              aVal with respect to the kth element of xState.
%        aHess A 7X7X7 collection of second partial derivatives of aVal
%              with respect to the elements of xState. ahess(:,k1,k2) is
%              the second partial derivative with respect to xState(k1) and
%              xState(k2).
%         papt The 7X1 vector of the partial derivative of aVal with
%              respect to time. This is all zeros.
%
%A derivation of the continuous-time flat-Earth circular turn model is
%given in [1] along with a description of how to use the model on a curved
%Earth.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(uTurn))
        uTurn=[0;0;1]; 
    end

    Omega=x(7)*uTurn;
    %We are evaluating
    % aDeriv=[x(4:6);cross(Omega,x(4:6));0];
    %But it is faster to do the cross product out element by element:
    aDeriv=[x(4:6);Omega(2)*x(6)-Omega(3)*x(5);Omega(3)*x(4)-Omega(1)*x(6);Omega(1)*x(5)-Omega(2)*x(4);0];
    
    if(nargout>1)
        dxdvy=-uTurn(3)*x(7);
        dxdvz=uTurn(2)*x(7);
        dxdx7=-uTurn(3)*x(5)+uTurn(2)*x(6);
        
        dydvx=uTurn(3)*x(7);
        dydvz=-uTurn(1)*x(7);
        dydx7=uTurn(3)*x(4)-uTurn(1)*x(6);
        
        dzdvx=-uTurn(2)*x(7);
        dzdvy=uTurn(1)*x(7);
        dzdx7=-uTurn(2)*x(4)+uTurn(1)*x(5);

        aJacob=[0,  0,    0,  1,        0,      0,      0;
                0,  0,    0,  0,        1,      0,      0;
                0,  0,    0,  0,        0,      1,      0;
                0,  0,    0,  0,        dxdvy,  dxdvz,  dxdx7;
                0,  0,    0,  dydvx,    0,      dydvz,  dydx7;
                0,  0,    0,  dzdvx,    dzdvy,  0,      dzdx7;
                0,  0,    0,  0,        0,      0,      0];

        if(nargout>2)
            aHess=zeros(7,7,7);
 
            dxdvydx7=-uTurn(3);
            dxdvzdx7=uTurn(2);
            dydvxdx7=uTurn(3);
            dydvzdx7=-uTurn(1);
            dzdvxdx7=-uTurn(2);
            dzdvydx7=uTurn(1);

            aHess(:,:,4)=[0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      dydvxdx7;
                          0,  0,    0,  0,        0,      0,      dzdvxdx7;
                          0,  0,    0,  0,        0,      0,      0];
            aHess(:,:,5)=[0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      dxdvydx7;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      dzdvydx7;
                          0,  0,    0,  0,        0,      0,      0];
            aHess(:,:,6)=[0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      dxdvzdx7;
                          0,  0,    0,  0,        0,      0,      dydvzdx7;
                          0,  0,    0,  0,        0,      0,      0;
                          0,  0,    0,  0,        0,      0,      0];
            aHess(:,:,7)=[0,  0,    0,  0,        0,                0,  0;
                          0,  0,    0,  0,        0,                0,  0;
                          0,  0,    0,  0,        0,                0,  0;
                          0,  0,    0,  0,        dxdvydx7,  dxdvzdx7,  0;
                          0,  0,    0,  dydvxdx7, 0,         dydvzdx7,  0;
                          0,  0,    0,  dzdvxdx7, dzdvydx7,         0,  0;
                          0,  0,    0,  0,        0,                0,  0];

            if(nargout>3)
                papt=zeros(7,1);
            end
        end
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
