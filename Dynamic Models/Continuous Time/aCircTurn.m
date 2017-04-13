function val=aCircTurn(x,t,uTurn)
%%ACIRCTURN The drift function for a constant-speed circular turning model
%           about an arbitrary axis.
%
%INPUTS: x The 7X1 target state at time t. It consists of 3D position,
%          velocity, and the scalar turn rate. The axis of the turn is a
%          parameter to this drift function.
%        t An unused time component so that aCircTurn can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
%    uTurn A unit vector indicating the axis of rotation about which the
%          turn is performed. If this parameter is omitted or an empty
%          matrix is passed, then uTurn=[0;0;1] is used, which means that
%          the turn takes place in the x-y plane.
%
%OUTPUTS: val The flat-Earth time-derivative of the state.
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

    if(nargin<3||isempty(uTurn))
        uTurn=[0;0;1]; 
    end

    Omega=x(7)*uTurn;
    val=[x(4:6);cross(Omega,x(4:6));0];
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
