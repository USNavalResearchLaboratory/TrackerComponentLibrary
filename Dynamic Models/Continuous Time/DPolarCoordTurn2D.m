function [D,DJacob,DHess,pDpt]=DPolarCoordTurn2D(x,qTurn,qSLin)
%%DPOLARCOORDTURN2D The continuous-time diffusion matrix function for a 2D
%              coordinated turn model where the velocity is specified as a
%              heading and an angle. The turn rate can be specified in
%              terms of a turn rate in radians per econd, or in terms of a
%              transversal acceleration. Additionally, a linear
%              acceleration can be given. This diffusion matrix goes with
%              the drift function aPolarCoordTurn2DOmega and
%              aPolarCoordTurn2DTrans.
%
%INPUTS: x The target state for 2D motion where the velocity is given in
%          terms of heading and speed components. If there is no linear
%          acceleration (acceleration along the direction of motion),
%          then x can either be x=[x;y;h;v;omega], where h is the heading
%          in terms of radians counterclockwise from the x-axis, v is the
%          speed, and omega is the turn rate (the derivative of h with
%          respect to time) or  x=[x;y;h;v;at] where at is the
%          transversal acceleration, which is orthogonal to
%          the velocity and is defined such that positive values of at
%          map to positive values of omega. If there is a linear
%          acceleration, then the target state is either
%          x=[x;y;h;v;omega;al] where omega is the turn rate and al
%          is the linear acceleration or the target state is
%          x=[x;y;h;v;at;al] if the turn is expressed in terms of a
%          transversal acceleration. The dimensionality of the state is
%          used to determine whether a linear acceleration component is
%          present. The linear acceleration component changes the speed.
%          That means that it is the derivative of the speed. The number of
%          columns in x determine how many copies of D are returned.
%    qTurn If the turn is specified in terms of a turn rate in radians
%          per second, then this is the power spectral density of the
%          turn rate noise having units of radians squared per seconds
%          cubed. If the turn is expressed in terms of a transverse
%          acceleration, then this is the power spectral density of the
%          transverse acceleration noise, having units of m^2/s^5.
%    qSLin If a linear acceleration is used, then this is the power
%          spectral density of the linear acceleration noise having units
%          of m^2/s^5. It a linear acceleration is not used, then this is
%          the power spectral density of the velocity noise having units of
%          m^2/s^3.
%
%OUTPUTS: D The xDimX2 diffusion matrix of a 2D continuous-time turning
%           model where the velocity is given as a heading and an angle. If
%           x has N columns, then N copies of D are returned with D(:,:,i)
%           being the ith one.
%   DJacob, DHess The xDimX2XxDim, xDimX2XxDimXxDim matrices of first and
%           second partial derivatives of the elements of D with respect to
%           x. These are all zero, since D is a constant. it is the same
%           for all D and is not repeated N times.
%      pDpt The xDimX2 partial derivative of D with respect to time. This
%           is all zeros, because D is a constant.
%
%The basic 2D polar coordinated turn model is described in [1]. It is also
%mentioned in [2], though no differential equations are given and a more
%detailed reference cited therein is a hard-to-get dissertation in French.
%The use of transversal acceleration is discussed in more detail in  [3].
%Chapter 4.2.3 of [4] also provides an overview of turning models using
%polar coordinates.
%
%More information on turning modeling using trasverse and linear
%acceleration can be found in the comments to the functions
%aCoordTurn2DOmega and aCoordTurn2DTrans.
%
%A starting point for setting qTurn when a turn rate is given and for
%setting qSLin when a linear acceleration is not given is to use
%processNoiseSuggest with 'PolyKal-ROT' and order=1. A starting point for
%setting qTurn when using a linear acceleration and for qSLin when a linear
%acceleration is given is to use processNoiseSuggest with 'PolyKal-ROT' and
%order=2. Similarly, The order chosen in the suggestion function just
%depends on the number of derivatives of time present.
%
%The corresponding drift functions are given by aPolarCoordTurn2DOmega and
%aPolarCoordTurn2DTrans. The corresponding discrete-time functions are
%FPolarCoordTurn2D and QPolarCoordTurn2D. However, note that the
%discrete-time functions with unknown noise is a direct-discrete model and
%not a discretization of the continuous-time model.
%
%REFERENCES:
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%[2] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[3] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design
%    of a multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%[4] S. Blackman and R. Popoli, Design and Analysis of Modern Tracking
%    Systems. Norwood, MA: Artech House, 1999.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(x,2);

rootQTurn=sqrt(qTurn);
rootQSLin=sqrt(qSLin);

numDim=length(x);
switch(numDim)
    case 5%If there is no linear acceleration
        D=[0,     0;%Row for x-component noise.
           0,     0;%Row for y-component noise.
           0,     0;%Row for heading noise.
           rootQSLin,0;%Row for speed noise.
           0,     rootQTurn];%Row for turn noise.
    case 6%If there is a linear acceleration term.
        D=[0,     0;%Row for x-component noise.
           0,     0;%Row for y-component noise.
           0,     0;%Row for heading noise.
           0,     0;%Row for speed noise.
           rootQTurn,0;%Row for turn noise.
           0,   rootQSLin];%Row for linear accelertion noise.
    otherwise
        error('The length of x is neither 5 nor 6.');
end

if(nargout>1)
    DJacob=zeros(numDim,2,numDim);
    if(nargout>2) 
        DHess=zeros(numDim,2,numDim,numDim);
        if(nargout>3)
            pDpt=zeros(numDim,2);
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
