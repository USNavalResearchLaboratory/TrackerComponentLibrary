function [D,DJacob,DHess,pDpt]=DCoordTurn2D(x,q0,qTurn,qLin)
%DCOORDTURN2D The continuous-time diffusion matrix function for a 2D
%             coordinated turn model with a Cartesian state. The turn rate
%             can be specified in terms of a turn rate in radians per
%             second, or in terms of a transversal acceleration.
%             Additionally, a linear acceleration can be given. This
%             diffusion matrix goes with the drift functions
%             aCoordTurn2DOmega and aCoordTurn2DTrans.
%
%INPUTS: x The target state for 2D motion. If there is no linear
%          acceleration (acceleration along the direction of motion), then
%          x can either be x=[x;y;xdot;ydot;omega], where omega is the turn
%          rate estimate in radians per second counterclockwise from the
%          x-axis or x=[x;y;xdot;ydot;at] where at is the transversal
%          acceleration, which is orthogonal to the velocity and is defined
%          such that positive values of at map to positive values of omega.
%          If there is a linear acceleration, then the target state is
%          either x=[x;y;xdot;ydot;omega;al] where omega is the turn rate
%          and al is the linear acceleration or the target state is
%          x=[x;y;xdot;ydot;at;al] if the turn is expressed in terms of a
%          transversal acceleration. The dimensionality of the state is
%          used to determine whether a linear acceleration component is
%          present. The linear acceleration component changes the speed.
%          That means that it acts in the direction of the velocity vector.
%          The number of columns in x determine how many copies of D are
%          returned.
%       q0 The power spectral density of the process noise of the velocity
%          components. It is assumed to be the same in both dimensions. It
%          covers perturbations from an ideal coordinated turn trajectory
%          and has units of m^2/s^3. If no process noise is desired for the
%          velocity (i.e. all perturbations should be covered by noise on
%          the turn component and on the linear acceleration), then a value
%          of zero should be passed.
%    qTurn If the turn is specified in terms of a turn rate in radians per
%          second, then this is the power spectral density of the turn rate
%          noise having units of radians squared per seconds cubed. If the
%          turn is expressed in terms of a transverse acceleration, then
%          this is the power spectral density of the transverse
%          acceleration noise, having units of m^2/s^5.
%     qLin This parameter is only needed if a linear acceleration is
%          present. It is the power spectral density of the linear
%          acceleration noise having units of m^2/s^5.
%
%OUTPUTS: D The diffusion matrix of a 2D continuous-time turning model
%           where the velocity is given as a Cartesian vector. If x has
%           N columns, then N copies of D are returned with D(:,:,i) being
%           the ith one.
%   DJacob, DHess The xDimXmXxDim, xDimXmXxDimXxDim matrices of first and
%           second partial derivatives of the elements of D with respect to
%           x. These are all zero, since D is a constant. it is the same
%           for all D and is not repeated N times.
%      pDpt The xDimX2 partial derivative of D with respect to time. This
%           is all zeros, because D is a constant.
%
%The basic 2D coordinated turn model in Cartesian coordinates is described
%in Section VA of [1]. When the turn rate is something that must be
%estimated, it is assumed that the continuous-time turn rate model is
%omegaDot=-(1/tau)*Omega+noise
%Note that the ordering of the state elements assumed by this function
%differs from the ordering of the state elements assumed in the paper.
%
%The 2D coordinates turn model in Cartesian coordinates is also described
%in  Chapter 4.2.3 of [2].
%
%The concept of using the transversal acceleration instead of the turn rate
%is not discussed in either of those references. It is, however, mentioned
%in [3], though no differential equations are given and a more detailed
%reference cited therein is a hard-to-get dissertation in French. The use
%of transversal acceleration is discussed in more detail in [4], though
%expressions are given when considering the 2D velocity are broken into
%components of heading and speed rather than in Cartesian space. The
%generalization to Cartesian space is not difficult and is done here.
%
%A starting point for setting q0 is to use processNoiseSuggest with
%'PolyKal-ROT' and order=1. The noise must cover small velocity deviations
%from the ideal turn model. A starting point for setting qTurn when a turn
%rate is given is to use processNoiseSuggest with 'PolyKal-ROT' and
%order=1. A starting point for setting qTurn when using a linear
%acceleration and for qLin is to use processNoiseSuggest with 'PolyKal-ROT'
%and order=2. The order chosen in the suggestion function just depends on
%the number of derivatives of time present.
%
%The corresponding drift functions are given by the functions
%aCoordTurn2DOmega and aCoordTurn2DTrans. The corresponding discrete-time
%functions are FCoordTurn2D and QCoordTurn. However, note that the
%discrete-time functions with unknown noise is a direct-discrete model and
%not a discretization of the continuous-time model.
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
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(x,2);

rootQ0=sqrt(q0);
rootQTurn=sqrt(qTurn);

numDim=length(x);
switch(numDim)
    case 5%There is no linear acceleration
        %Equation 61 in Li's paper combined with Equation 67.
        D=[0,       0,      0;%Row for x-component noise.
           0,       0,      0;%Row for y-component noise.
           rootQ0,  0,      0;%Row for velocity-x noise.
           0,       rootQ0, 0;%Row for velocity-y noise.
           0,       0,      rootQTurn];%Row for turn noise.
       m=3;
    case 6%There is a linear acceleration component.
        rootQLin=sqrt(qLin);
        %Similar to Equation 61 in Li's paper combined with Equation 67,
        %but with an added row for noise in the linear acceleration
        %component.
        D=[0,       0,      0,  0;%Row for x-component noise.
           0,       0,      0,  0;%Row for y-component noise.
           rootQ0,  0,      0,  0;%Row for velocity-x noise.
           0,       rootQ0, 0,  0;%Row for velocity-y noise.
           0,       0, rootQTurn,0;%Row for turn rate noise.
           0,       0,      0,  rootQLin];%Row for linear accel. noise.
       m=4;
    otherwise
        error('The length of x is neither 5 nor 6.');
end

if(nargout>1)
    DJacob=zeros(numDim,m,numDim);
    if(nargout>2) 
        DHess=zeros(numDim,m,numDim,numDim);
        if(nargout>3)
            pDpt=zeros(numDim,m);
        end
    end
end

if(N>1)
   D=repmat(D,[1,1,N]); 
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
