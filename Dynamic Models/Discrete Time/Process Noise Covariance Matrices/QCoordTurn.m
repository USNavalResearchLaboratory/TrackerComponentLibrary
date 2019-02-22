function Q=QCoordTurn(T,x,sigmaV2,sigmaTurn2,sigmaLin2)
%%QCOORDTURN Get the discrete-time process noise covariance matrix for a
%            2D or 3D coordinated turn model with a Cartesian state. The
%            turn rate can be specified in terms of a turn rate in radians
%            per second, or in terms of a transversal acceleration.
%            Additionally, a linear acceleration can be given.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%        x The target state for 2D or 3D motion. If there is no linear
%          acceleration (acceleration along the direction of motion), then
%          x in 2D can either be x=[x;y;xdot;ydot;omega], where omega is
%          the turn rate estimate in radians per second counterclockwise
%          from the x-axis or x can be x=[x;y;xdot;ydot;at] where at is
%          the transversal acceleration, which is orthogonal to the
%          velocity and is defined such that positive values of at map to
%          positive values of omega. In 3D, the state has z and zdot
%          components after y and ydot. If there is a linear acceleration,
%          then in 2D the target state is either x=[x;y;xdot;ydot;omega;al]
%          where omega is the turn rate and al is the linear acceleration
%          or the target state is x=[x;y;xdot;ydot;at;al] if the turn is
%          expressed in terms of a transversal acceleration. Again, in 3D,
%          z and zdot components are respectively added after y and ydot. The
%          dimensionality of the state is used to determine the
%          dimensionality and whether a linear acceleration component is
%          present. The linear acceleration component changes the speed.
%          That means that it acts in the direction of the velocity vector.
%  sigmaV2 The variance driving the process noise affecting the velocity
%          components. This has units of m^2/s^4 and is assumed to be the
%          same in both the x and y dimensions.
% sigmaTurn2 If the turn is specified in terms of a turn rate in radians
%          per second, then this is the variance driving the process noise
%          of the turn rate having units of radians squared per second
%          squared. If the turn is expressed in terms of a transverse
%          acceleration, then this is the variance of the transverse
%          acceleration noise, having units of m^2/s^4.
% sigmaLin2 This parameter is only needed if a linear acceleration is
%          present. It is the variance of the linear acceleration noise
%          having units of m^2/s^4.
%
%OUTPUTS: Q The process noise covariance matrix under a direct discrete-
%           time coordinated turn model where the velocity is specified in
%           Cartesian coordinates. Note that Q is singular, due to the
%           direct-discrete modeling of the process noise.
%
%The basic 2D coordinated turn model is described in Section VA of [1]. It
%is assumed that the continuous-time turn rate model is
%omegaDot=-(1/tau)*Omega+noise, which discretizes to
%omega[k+1]=exp(-T/tau)*omega[k]+noise.
%The turn rate for an unknown process noise is taken from Equation 73 of 
%the paper and is a discrete-time approximation, not an exact
%discretization. The ordering of the elements in the state has been changed
%from the paper.
%
%The velocity and position components of the process noise model are
%similar to a discrete white noise acceleration model. Thus, a starting
%point for the process noise parameter sigmaV2 can be found using the
%processNoiseSuggest function with 'PolyKalDirectDisc-ROT' as the model and
%order=1 to cover velocity changes that are not solely due to turns or
%linear accelerations. The changes in the turn rate/transverse acceleration
%that are not due to exponential decay can similarly be covered using
%process noise  starting with sigmaTurn2 being set based on
%processNoiseSuggest with 'PolyKalDirectDisc-ROT' and order=1/2. Similarly,
%a starting point for sigmaLin2 is from processNoiseSuggest with
%'PolyKalDirectDisc-ROT' and order=2.
%
%Note that the process noise added to the turn component and the linear
%acceleration term (if present) is not integrated into the velocity/
%position components. If the only changes allowed/ expected are due to
%modifications in the linear and transverse accelerations (or turn rate),
%then sigmaV2=0. However, regardless of the size of T, the uncertainty in
%the other components would not directly integrate into the position and
%velocity components over any one step. Thus, one would expect the model to
%be particularly poor when the step size is large.
%
%The corresponding transition matrix in 2D is given by the function
%FCoordTurn2D. The corresponding continuous-time functions are
%aCoordTurn2DOmega and aCoordTurn2DTrans with diffusion matrix
%DCoordTurn2D. However, note that the discrete-time functions with
%unknown noise is a direct-discrete model and not a discretization of the
%continuous-time model.
%
%REFERENCES:
%[1] X. R. Li and V. P. Jilkov, "Survey of maneuvering target tracking.
%   Part I: Dynamic models," IEEE Transactions on Aerospace and Electronic
%   Systems, vol. 39, no. 4, pp. 1333-1364, Oct. 2003.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

sigmaV=sqrt(sigmaV2);
sigmaTurn=sqrt(sigmaTurn2);

switch(length(x))
    case 5%2D, there is no linear acceleration.        
        G=[T^2/2;%Position x
           T^2/2;%Position y
           T;%Velocity x
           T;%Velocity y
           T];%Turn rate/ linear acceleration
        G(1:4)=G(1:4)*sigmaV;
        G(5)=G(5)*sigmaTurn;
        Q=G*G';
    case 6%2D, there is linear acceleration.
        sigmaLin=sqrt(sigmaLin2);
        G=[T^2/2;%Position x
           T^2/2;%Position y
           T;%Velocity x
           T;%Velocity y
           T;%Turn rate/ linear acceleration
           T];%Transverse acceleration
        G(1:4)=G(1:4)*sigmaV;
        G(5)=G(5)*sigmaTurn;
        G(6)=G(6)*sigmaLin;
        Q=G*G';
    case 7%3D, there is no linear acceleration.        
        G=[T^2/2;%Position x
           T^2/2;%Position y
           T^2/2;%Position z
           T;%Velocity x
           T;%Velocity y
           T;%Velocity z
           T];%Turn rate/ linear acceleration
        G(1:6)=G(1:6)*sigmaV;
        G(7)=G(7)*sigmaTurn;
        Q=G*G';
    case 8%3D, there is linear acceleration.
        sigmaLin=sqrt(sigmaLin2);
        G=[T^2/2;%Position x
           T^2/2;%Position y
           T^2/2;%Position z
           T;%Velocity x
           T;%Velocity y
           T;%Velocity z
           T;%Turn rate/ linear acceleration
           T];%Transverse acceleration
        G(1:6)=G(1:6)*sigmaV;
        G(7)=G(7)*sigmaTurn;
        G(8)=G(8)*sigmaLin;
        Q=G*G';
    otherwise
        error('The length of x is neither 5 nor 6 (for 2D) and neither 7 or 8 (for 3D).');
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
