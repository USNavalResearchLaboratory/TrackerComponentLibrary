function val=DPolarLin2D(x,t,q0Dir,param4)
%%DPOLARCV2D The diffusion matrix for a continuous-time constant heading
%            motion model where the target state is given in terms of
%            position, direction of velocity in radians counterclockwise
%            from the x-axis, and speed.
%
%INPUTS: x The 4XN set of target state vectors in 2D space in the order of
%          [2D position;direction angle;speed] for each column where the
%          heading angle is measured in radians counterclockwise from the
%          x-axis. Alternatively, a state vector of
%          [2D position;direction angle;speed; speed derivative] can be
%          used to cover a linearlly accelerating target. Only the
%          dimensionality of x matters; the components are ignored. The
%          dimensionality is used to determine whether linear acceleration
%          is modeled and what N is.
%        t An unused time component so that aLinear can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
%    q0Dir The power spectral density of the process noise for the heading
%          direction. The units are rad^2/s.
%  param4  If x is only 4D, that is, does not model a constant linear
%          acceleration, then param4 is the power spectral density of the
%          process noise for the speed. The units are m^2/s^3. Otherwise,
%          param4 is the power spectral density of the noise for the
%          linear acceleration and its units are m^2/s^5.
%
%OUTPUT: D    The 4X2XN set of N diffusion matrices of a continuous-time
%             linear additive noise model where the velocity is given in
%             terms of direction and speed and the noise is only added to
%             the direction and speed components or the 5X2XN set of
%             matrices whrere noise is added to the direction and
%             acceleration components if a linear acceleration is assumed.
%
%Ideas for setting q0Dir  can be obtained using processNoiseSuggest with
%'PolyKal-ROT' and order 0. Ideas for setting param4 can
%be obtained using processNoiseSuggest with 'PolyKal-ROT' and order 1 or 2
%respectively depending on whether x is constant velocity or constant
%acceleration.
%
%The idea of decomposing the velocity into direction and speed components
%is presented in [1], where the direction component is differently defined
%and additional derivatives are used. 
%
%This diffusion matrix goes with the drift function aPolarLin2D.
%
%REFERENCES:
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(x,2);

q0DirR=sqrt(q0Dir);
param4R=sqrt(param4);

if(size(x,1)==4)%If it is a constant velocity model
    val=repmat([0,      0;
                0,      0;
                q0DirR, 0;
                0,     param4R],[1,1,N]);
else%If it is a constant acceleration model.
    val=repmat([0,      0;
                0,      0;
                q0DirR, 0;
                0,      0;
                0,      param4R],[1,1,N]);
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
