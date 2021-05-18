function xPol=Cart2DState2PolarState(xCart,systemType)
%%CART2DSSTATE2POLARSTATE Transform a 2D Cartesian state into a state
%                         consisting of position, heading and speed as well
%                         as possibly a turn rate and a linear
%                         acceleration, depending on the choice of
%                         systemType.
%
%INPUTS: xCart A Cartesian state vector consisting of position velocity and
%              possibly acceleration into a state where heading and speed
%              have been separated. xCart has the form
%              [x;y;xdot;ydot;xddot;yddot], where the acceleration terms
%              xddot;yddot can be omitted if the system type is 'ConstVel'.
%   systemType A string constant specifying the desired type of output. In
%              all instances, the heading is measured in terms of radians
%              counterclockwise from the x-axis. Possible values are:
%              'ConstVel'     The target state is [position;heading;speed]
%                             and xCart is [position;velocity]
%              'ConstAccel'   The target state is [position;heading;speed;
%                             speed derivative] and xCart is
%                             [position;velocity;acceleration]
%              'ConstTurn'    The target state is [position;heading;speed;
%                             turn rate] and xCart is
%                             [position;velocity;acceleration]
%              'TurnAndAccel' The target state is [position;heading;speed;
%                             turnrate; speed derivative] and xCart is
%                             [position;velocity;acceleration]
%
%OUTPUTS: xPol The state converted from 2D Cartesian coordinates into the
%              selected 2D coordinate system.
%
%When the system type is 'ConstVel' or 'TurnAndAccel', only a single
%solution is mathematically observable. When the system type is
%'ConstAccel' or 'ConstTurn', the system is overdetermined, but only a
%simple solution is used, not a least squares solution.
%
%The use of 2D states where the heading and speed have been separated is
%discussed in [1] and [2].
%
%The opposite of this function is polar2DState2CartState.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get position
x=xCart(1);
y=xCart(2);
%Get velocity
xDot=xCart(3);
yDot=xCart(4);

switch(systemType)
    case 'ConstVel' %Given position, heading, and speed.
        xPol=[x;
              y;
              atan2(yDot,xDot);
              sqrt(yDot^2+xDot^2)];
    case 'ConstAccel'%Given position, heading, speed, and linear
                     %acceleration
        %Get acceleration
        xDdot=xCart(5);
        yDdot=xCart(6);
        
        theta=atan2(yDot,xDot);%Heading
        %Determine the sign of the derivative velocity, ignoring the possible
        %effects of noise...
        vDot=sqrt(yDdot^2+xDdot^2);%Linear acceleration
        
        diff1=(vDot*cos(theta)-xDdot)^2+(vDot*sin(theta)-yDdot)^2;
        diff2=(-vDot*cos(theta)-xDdot)^2+(-vDot*sin(theta)-yDdot)^2;
        if(diff2<diff1)
            vDot=-vDot;
        end
        xPol=[x;
              y;
              theta;
              sqrt(yDot^2+xDot^2);
              vDot];
    case 'ConstTurn' %Given position, heading, speed, and turn rate.
        %Get acceleration
        xDdot=xCart(5);
        yDdot=xCart(6);
    
        %Turn rate
        omega=(xDot*yDdot-yDot*xDdot)/(xDot^2+yDot^2);

        if(~isfinite(omega))
            omega=0;
        end

        xPol=[x;
              y;
              atan2(yDot,xDot);
              sqrt(yDot^2+xDot^2);
              omega];
    case 'TurnAndAccel'%Given position, heading, speed, turn rate, and
                       %linear acceleration.
        %Get acceleration
        xDdot=xCart(5);
        yDdot=xCart(6);
        
        theta=atan2(yDot,xDot);%Heading
        v=sqrt(yDot^2+xDot^2);%Speed
        omega=(yDdot*cos(theta)-xDdot*sin(theta))/v;%Turn rate

        %Deal with slow speed targets. 
        if(~isfinite(omega))
            omega=0;
        end

        %Linear acceleration
        vDot=xDdot*cos(theta)+yDdot*sin(theta);

        xPol=[x;
              y;
              theta;
              v;
              omega;
              vDot];
    otherwise
        error('Invalid system type given.')
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
