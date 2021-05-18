function [aVal,aJacob,aHess,papt]=aCoordTurn2DOmega(xState,tauInv,tauAccelInv)
%%ACCORDTURN2DOMEGA The continuous-time drift function for a 2D coordinated
%              turn model with a Cartesian state and a turn rate given as
%              part of the state. Additionally, a linear acceleration term
%              can be given. It acts in the direction of motion. The turn
%              rate and linear acceleration can optionally have time
%              constants associated with them, like in the Singer model,
%              modelling a tendancy to eventually return to non-
%              accelerating, straight-line motion.
%
%INPUTS: xState The 5X1 or 6X1 target state for 2D motion. If there is no
%               linear acceleration (acceleration along the direction of
%               motion), then xState=[x;y;xdot;ydot;omega], where omega is
%               the turn rate counterclockwise from the x axis (typically
%               radians per second).  If there is a linear acceleration,
%               then the target state is xState=[x;y;xdot;ydot;omega;al],
%               where al is the linear acceleration. The dimensionality of
%               the state is used to determine whether a linear
%               acceleration component is present. The linear acceleration
%               component changes the speed. That means that it acts in the
%               direction of the velocity vector.
%       tauInv The inverse of the correlation time constant tau for the
%              turn rate in seconds. tauInv must be positive. The default
%              if omitted or an empty matrix is passed is zero.
%  tauAccelInv The inverse of the correlation time constant for the linear
%              acceleration in seconds. This parameter is not used if there
%              is no linear acceleration. the default if omitted or an
%              empty matrix is passed is 0.
%
%OUTPUTS: aVal The 5X1 (or 6X1 with linear acceleration) time-derivative of
%              the state. 
%       aJacob The 5X5 (or 6X6) matrix of partial derivatives of aVals
%              such that aJacob(:,k) is the partial derivative of
%              aVals(:,k) with respect to xState(k).
%        aHess The 5X5X5 (or 6X6X6) matrix of second derivatives of aVals
%              such that aHess(:,k1,k2) is the second partial derivative of
%              aVals with respect to xState(k1) and xState(k2).
%         papt The 5X1 or 6X1 derivative with resect to time of aVals.
%              This is all zeros, because the model is time invariant.
%
%The basic 2D coordinated turn model in Cartesian coordinates is described
%in Section VA of [1]. When the turn rate is something that must be
%estimated, it is assumed that the continuous-time turn rate model is
%omegaDot=-(1/tauTurn)*omega+noise
%Note that the ordering of the state elements assumed by this function
%differs from the ordering of the state elements assumed in [1].
%
%The 2D coordinates turn model in Cartesian coordinates is also described
%in Chapter 4.2.3 of [2].
%
%The concept of using the transversal acceleration instead of the turn rate
%is not discussed in either of those references. It is, however, mentioned
%in [3], though no differential equations are given and a more detailed
%reference cited therein is a hard-to-get dissertation in French. The use
%of transversal acceleration is discussed in more detail in [4], though
%expressions are given when considering the 2D velocity are broken
%into components of heading and speed rather than in Cartesian space. The
%generalization to Cartesian space is not difficult and is done here.
%
%The corresponding diffusion matrix is given by the function DCoordTurn2D.
%The corresponding discrete-time functions are FCoordTurn2D and
%QCoordTurn. However, note that the discrete-time functions are
%direct-discrete models and not discretizations of the continuous-time
%models as the propagated PDF does not remain Gaussian over time.
%
%REFERENCES:
%[1] X. R. Li and V. P. Jilkov, "Survey of maneuvering target tracking.
%    Part I: Dynamic models," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 39, no. 4, pp. 1333-1364, Oct. 2003.
%[2] S. Blackman and R. Popoli, Design and Analysis of Modern Tracking
%    Systems. Norwood, MA: Artech House, 1999.
%[3] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[4] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design
%    of a multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(tauAccelInv))
    tauAccelInv=0; 
end

if(nargin<2||isempty(tauInv))
    tauInv=0; 
end

xDot=xState(3);
yDot=xState(4);
omega=xState(5);

switch(length(xState))
    case 5%There is no linear acceleration. This model is from Equations 61
           %and 67 in [1]
        aVal=[xDot;%Position dervative
              yDot;%Position derivative
              -omega*yDot;%Velocity derivative
              omega*xDot;%Velocity derivative
              -tauInv*omega];%Turn rate derivative

        if(nargout>1)
            aJacob=[0,0,1,      0,      0;
                    0,0,0,      1,      0;
                    0,0,0,      -omega, -yDot;
                    0,0,omega,  0,      xDot;
                    0,0,0,      0,      -tauInv];
            if(nargout>2)
                aHess=zeros(5,5,5);
                aHess(:,:,3)=[0,0,0,0,0;
                              0,0,0,0,0;
                              0,0,0,0,0;
                              0,0,0,0,1;
                              0,0,0,0,0]; 
                aHess(:,:,4)=[0,0,0,0, 0;
                              0,0,0,0, 0;
                              0,0,0,0,-1;
                              0,0,0,0, 0;
                              0,0,0,0, 0]; 
                aHess(:,:,5)=[0,0,0,0, 0;
                              0,0,0,0, 0;
                              0,0,0,-1,0;
                              0,0,1,0, 0;
                              0,0,0,0, 0];
                if(nargout>3)
                    papt=zeros(5,1); 
                end
            end
        end
    case 6
        al=xState(6);
        
        xDot2=xDot*xDot;
        yDot2=yDot*yDot;
        
        v2=xDot2+yDot2;
        v=sqrt(v2);%The speed

        %Get the linear acceleration vector. There will be NaNs if
        %the velocity is zero.
        linAccel=([xDot;yDot]/v)*al;

        %Deal with problems with a zero or nearly zero velocity.
        if(any(~isfinite(linAccel)))
            zeroTerms=true;
            linAccel(:)=0;
        else
            zeroTerms=false;
        end
        
        %From Equations 61 and 67 in Li's paper with linear
        %acceleration added.
        aVal=[xDot;%Position dervative
              yDot;%Position derivative
              -omega*yDot+linAccel(1);%Velocity derivative
              omega*xDot+linAccel(2);%Velocity derivative
              -tauInv*omega;%Turn rate derivative
              -tauAccelInv*al];%Linear acceleration derivative

        if(nargout>1)
            if(zeroTerms==false)
                v3=v*v2;
                plaxpxDot=(al*yDot2)/v3;
                plaxpyDot=-((al*xDot*yDot)/v3);
                plaxpal=xDot/v;

                playpxDot=-((al*xDot*yDot)/v3);
                playpyDot=(al*xDot2)/v3;
                playpal=yDot/v;
            else
                plaxpxDot=0;
                plaxpyDot=0;
                plaxpal=0;
                playpxDot=0;
                playpyDot=0;
                playpal=0;
            end

            aJacob=[0,0,1,              0,                  0,      0;
                    0,0,0,              1,                  0,      0;
                    0,0,plaxpxDot,      -omega+plaxpyDot,   -yDot,  plaxpal;
                    0,0,omega+playpxDot,playpyDot,          xDot,   playpal;
                    0,0,0,              0,                  -tauInv,0;
                    0,0,0,              0,                  0,      -tauAccelInv];

            if(nargout>2)
                if(zeroTerms==false)
                    v5=v3*v2;
                    plaxpxDotpxDot=-((3*al*xDot*yDot2)/v5);
                    plaxpxDotpyDot=-((al*yDot*(-2*xDot2+yDot2))/v5);
                    plaxpxDotpal=yDot^2/v3;

                    plaxpyDotpxDot=plaxpxDotpyDot;
                    plaxpyDotpyDot=-((al*xDot*(xDot2-2*yDot2))/v5);
                    plaxpyDotpal=-((xDot*yDot)/v3);

                    plaxpalpxDot=plaxpxDotpal;
                    plaxpalpyDot=plaxpyDotpal;

                    playpxDotpxDot=-((al*yDot*(-2*xDot2+yDot2))/v5);
                    playpxDotpyDot=-((al*xDot*(xDot2-2*yDot2))/v5);
                    playpxDotpal=-((xDot*yDot)/v3);

                    playpyDotpxDot=playpxDotpyDot;
                    playpyDotpyDot=-((3*al*xDot2*yDot)/v5);
                    playpyDotpal=xDot^2/v3;

                    playpalpxDot=playpxDotpal;
                    playpalpyDot=playpyDotpal;
                else
                    plaxpxDotpxDot=0;
                    plaxpxDotpyDot=0;
                    plaxpxDotpal=0;
                    plaxpyDotpxDot=0;
                    plaxpyDotpyDot=0;
                    plaxpyDotpal=0;
                    plaxpalpxDot=0;
                    plaxpalpyDot=0;
                    playpxDotpxDot=0;
                    playpxDotpyDot=0;
                    playpxDotpal=0;
                    playpyDotpxDot=0;
                    playpyDotpyDot=0;
                    playpyDotpal=0;
                    playpalpxDot=0;
                    playpalpyDot=0;
                end

                aHess=zeros(6,6,6);
                aHess(:,:,3)=[0,0,0,              0,                  0,      0;
                              0,0,0,              0,                  0,      0;
                              0,0,plaxpxDotpxDot, plaxpyDotpxDot,     0,      plaxpalpxDot;
                              0,0,playpxDotpxDot, playpyDotpxDot,     1,      playpalpxDot;
                              0,0,0,              0,                  0,      0;
                              0,0,0,              0,                  0,      0];
                aHess(:,:,4)=[0,0,0,              0,                  0,      0;
                              0,0,0,              0,                  0,      0;
                              0,0,plaxpxDotpyDot, plaxpyDotpyDot,     -1,     plaxpalpyDot;
                              0,0,playpxDotpyDot, playpyDotpyDot,     0,      playpalpyDot;
                              0,0,0,              0,                  0,      0;
                              0,0,0,              0,                  0,      0];
                aHess(:,:,5)=[0,0,0, 0,0,0;
                              0,0,0, 0,0,0;
                              0,0,0,-1,0,0;
                              0,0,1, 0,0,0;
                              0,0,0, 0,0,0;
                              0,0,0, 0,0,0];
                aHess(:,:,6)=[0,0,0,              0,              0,0;
                              0,0,0,              0,              0,0;
                              0,0,plaxpxDotpal,   plaxpyDotpal,   0,0;
                              0,0,playpxDotpal,   playpyDotpal,   0,0;
                              0,0,0,              0,              0,0;
                              0,0,0,              0,              0,0];
                if(nargout>3)
                    papt=zeros(6,1); 
                end
            end
        end
    otherwise
        error('The dimensionality of the state is neither 5 nor 6.');
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
