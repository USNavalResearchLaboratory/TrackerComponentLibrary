function val=aPolarLin2D(x,t,T)
%APOLARLIN2D The drift function for a 2D continuous-time constant heading
%            motion model where the target state is given in terms of
%            position, direction of velocity in radians counterclockwise
%            from the x-axis, and speed. Optionally, a speed derivative
%            component can be given to model possible linear target
%            acceleration.
%
%INPUTS: x The 4XN set of target state vectors in 2D space in the order of
%          of [2D position;direction angle;speed] where the heading angle
%          is measured in radians counterclockwise from the x-axis.
%          Alternatively, a 5XN state vector of [2D position;direction
%          angle;speed; speed derivative] can be used to cover a linearly
%          accelerating target.
%        t An unused time component so that aPolar2DCV can be used with
%          Runge-Kutta methods that expect the function to take two
%          parameters.
%        T The time-duration of the propagation interval in seconds. This
%          is only needed when linear acceleration is included.
%
%OUTPUTS: val The 4XN or 5XN set of time-derivatives of the state under the
%             constant heading motion model.
%
%The idea of decomposing the velocity into direction and speed components
%is presented in [1], where the direction component is differently defined
%and derivative for a turning model are used. 
%
%This drift function goes with the diffusion matrix DPolarLin2D.
%
%REFERENCES:
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(x,2);
val=zeros(size(x));

for curX=1:N
    %The direction component.
    theta=x(3,curX);
    v=x(4,curX);

    if(size(x,1)==4)%If it is just constant velocity.
        val(:,curX)=[v*cos(theta);
             v*sin(theta);
             0;
             0];
    else%If a linear acceleration is also given.
        vDot=x(5,curX);
        val(:,curX)=[(v+(T/2)*vDot)*cos(theta);
             (v+(T/2)*vDot)*sin(theta);
             0;
             vDot;
             0];
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
