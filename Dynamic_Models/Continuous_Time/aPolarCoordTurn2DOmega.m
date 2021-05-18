function [aVal,aJacob,aHess,papt]=aPolarCoordTurn2DOmega(x,tauInv,tauAccelInv)
%%APOLARCOORDTURN2DOMEGA The continuous-time drift function for a 2D
%           coordinated turn model with the velocity expressed in terms of
%           a heading angle and a speed rather than in terms of a Cartesian
%           velocity vector. Additionally a turn rate is part of the state
%           (unlike in aPolarLin2D) and a linear acceleration term can be
%           given, which acts in the direction of motion. The turn rate and
%           linear acceleration can optionally have time constants
%           associated with them, like in the Singer model, modelling a
%           tendancy to eventually return to non-accelerating, straight-
%           line motion.
%
%INPUTS: xState The 5X1 or 6X1 target state for 2D motion. If there is no
%               linear acceleration (acceleration along the direction of
%               motion), then the state is [xPos;yPos;theta;v;omega],
%               where (xPos,yPos) are the Cartesian position, theta is the
%               angle of the velocity vector in radians counterclockwise
%               from the x axis, v is the speed, and omega is the turn
%               rate, typically in radians per second. If a linear
%               acceleration component is provided, then the state is
%               [xPos;yPos;theta;v;omega,al].
%        tauInv The inverse of the correlation time constant tau for the
%               turn rate in seconds. tauInv must be positive. The default
%               if omitted or an empty matrix is passed is 0.
%   tauAccelInv The inverse of the correlation time constant for the linear
%               acceleration in seconds. This parameter is not used if
%               there is no linear acceleration. the default if omitted or
%               an empty matrix is passed is 0.
%
%OUTPUTS: aVal The 5X1 (or 6X1 with linear acceleration) time-derivative of
%              the state. 
%       aJacob The 5X5 (or 6X6) matrix of partial derivatives of aVals
%              such that aJacob(:,k) is the partial derivative of
%              aVals(:,k) with respect to x(k).
%        aHess The 5X5X5 (or 6X6X6) matrix of second derivatives of aVals
%              such that aHess(:,k1,k2) is the second partial derivative of
%              aVals with respect to x(k1) and x(k2).
%         papt The 5X1 or 6X1 derivative with resect to time of aVals.
%              This is all zeros, because the model is time invariant.
%
%The basic 2D coordinated turn model in Cartesian coordinates is described
%in Section VA of [1]. When the turn rate is something that must be
%estimated, it is assumed that the continuous-time turn rate model is
%omegaDot=-(1/tauTurn)*Omega+noise
%Note that the ordering of the state elements assumed by this function
%differs from the ordering of the state elements assumed in [1].
%
%The 2D coordinates turn model in Cartesian coordinates is also described
%in Chapter 4.2.3 of [2].
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
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(tauInv))
    tauInv=0;
end

if(nargin<3||isempty(tauAccelInv))
    tauAccelInv=0;
end

%Extract the heading and speed (polar velocity components).
theta=x(3);
v=x(4);
omega=x(5);%The turn rate

xDim=size(x,1);

sinTheta=sin(theta);
cosTheta=cos(theta);
switch(xDim)
    case 5%There is no linear acceleration
        %From Equations 61 and 67 in Li's paper.
        aVal=[v*cosTheta;%Position dervative
              v*sinTheta;%Position derivative
              omega;%Heading derivative
              0;%Speed derivative
              -tauInv*omega];%Turn rate derivative
        
        if(nargout>1)
            dXdTheta=-v*sinTheta;
            dYdTheta=v*cosTheta;
            dXdv=cosTheta;
            dYdv=sinTheta;
            
            aJacob=[0,0,dXdTheta,dXdv,0;
                    0,0,dYdTheta,dYdv,0;
                    0,0,0,          0,1;
                    0,0,0,          0,0;
                    0,0,0,          0,-tauInv];
            
            if(nargout>2)
                aHess=zeros(5,5,5);
            
                dXdThetadTheta=-v*cosTheta;
                dXdvdTheta=-sinTheta;
                
                dXdThetadv=dXdvdTheta;
                dXdvdv=0;
                
                %%%
                dYdThetadTheta=-v*sinTheta;
                dYdvdTheta=cosTheta;
                
                dYdThetadv=dYdvdTheta;
                dYdvdv=0;
                
                aHess(:,:,3)=[0,0,dXdThetadTheta,dXdvdTheta,0;
                              0,0,dYdThetadTheta,dYdvdTheta,0;
                              zeros(3,5)];
                aHess(:,:,4)=[0,0,dXdThetadv,dXdvdv,0;
                              0,0,dYdThetadv,dYdvdv,0;
                              zeros(3,5)];

                if(nargout>3)
                    papt=zeros(5,1);
                end
            end
        end
    case 6%There is a linear acceleration component.
        al=x(6);%The linear acceleration (vDot)
        
        aVal=[v*cosTheta;%Position dervative
              v*sinTheta;%Position derivative
              omega;%Heading derivative
              al;%Speed derivative
              -tauInv*omega;%Turn rate derivative
              -tauAccelInv*al];%Linear acceleration derivative
          
        if(nargout>1)
            dXdTheta=-v*sinTheta;
            dYdTheta=v*cosTheta;
            dXdv=cosTheta;
            dYdv=sinTheta;

            aJacob=[0,0,dXdTheta,dXdv,        0,0;
                    0,0,dYdTheta,dYdv,        0,0;
                    0,0,       0,   0,        1,0;
                    0,0,       0,   0,        0,1;
                    0,0,       0,   0,  -tauInv,0;
                    0,0,       0,   0,        0,-tauAccelInv];
            
            if(nargout>2)
                aHess=zeros(6,6,6);
            
                dXdThetadTheta=-v*cosTheta;
                dXdvdTheta=-sinTheta;
                
                dXdThetadv=dXdvdTheta;
                dXdvdv=0;
                
                %%%
                dYdThetadTheta=-v*sinTheta;
                dYdvdTheta=cosTheta;
                
                dYdThetadv=dYdvdTheta;
                dYdvdv=0;
                
                aHess(:,:,3)=[0,0,dXdThetadTheta,dXdvdTheta,0,0;
                              0,0,dYdThetadTheta,dYdvdTheta,0,0;
                              zeros(4,6)];

                aHess(:,:,4)=[0,0,dXdThetadv,dXdvdv,0,0;
                              0,0,dYdThetadv,dYdvdv,0,0;
                              zeros(4,6)];

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
