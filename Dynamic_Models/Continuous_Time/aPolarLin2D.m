function [aVal,aJacob,aHess,papt]=aPolarLin2D(x,T)
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
%        T The time-duration of the propagation interval in seconds. This
%          is only needed when linear acceleration is included.
%
%OUTPUTS: aVal The 4XN or 5XN set of time-derivatives of the state under
%              the constant heading motion model.
%       aJacob This and higher partial derivatives can only be requested
%              if N=1. This is the 4X4 or 5X5 matrix of partial derivatives
%              of aVal such that aJacob(:,i) is the partial derivative of
%              aVal with respect to x(i).
%        aHess The 4X4X4 or 5X5X5 matrix of second derivatives of aVal such
%              that aHess(:,k1,k2) is the second partial derivative of
%              aVal with respect to x(k1) and x(k2).
%         papt The 4X1 or 5X1 partial derivative with resect to time of
%              aVal. This is all zeros, because the model is time
%              invariant.
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
aVal=zeros(size(x));

for curX=1:N
    %The direction component.
    theta=x(3,curX);
    v=x(4,curX);

    cosTheta=cos(theta);
    sinTheta=sin(theta);
    
    if(size(x,1)==4)%If it is just constant velocity.
        aVal(:,curX)=[v*cosTheta;
             v*sinTheta;
             0;
             0];
         
        if(nargout>1)
            if(N>1)
                error('Derivatives are only available for N=1.')
            end
            aJacob=[0,0,-v*sinTheta,   cosTheta;
                    0,0,v*cosTheta,    sinTheta;
                    0,0,0,             0;
                    0,0,0,             0];
                
            if(nargout>2)
                aHess=zeros(4,4,4);
                
                aHess(:,:,3)=[0,0,-v*cosTheta,-sinTheta;
                              0,0,-v*sinTheta,cosTheta;
                              zeros(2,4)];
                aHess(:,:,4)=[0,0,-sinTheta,0;
                              0,0,cosTheta, 0;
                              zeros(2,4)];
                if(nargout>3)
                    papt=zeros(4,1);
                end     
            end
        end
    else%If a linear acceleration is also given.
        vDot=x(5,curX);
        
        vSum=(v+(T/2)*vDot);
        T2=T/2;
        
        aVal(:,curX)=[vSum*cosTheta;
                      vSum*sinTheta;
                      0;
                      vDot;
                      0];
        if(nargout>1)
            aJacob=[0,0,-vSum*sinTheta, cosTheta,T2*cosTheta;
                    0,0,vSum*cosTheta,  sinTheta,T2*sinTheta;
                    0,0,0,              0,        0;
                    0,0,0,              0,        1;
                    0,0,0,              0,        0];
            if(nargout>2)
                aHess=zeros(5,5,5);
                
                aHess(:,:,3)=[0,0,-vSum*cosTheta,-sinTheta,-T2*sinTheta;
                              0,0,-vSum*sinTheta,cosTheta,T2*cosTheta;
                              zeros(3,5)];
                aHess(:,:,4)=[0,0,-sinTheta,0,0;
                              0,0,cosTheta,0,0;
                              zeros(3,5)];
                aHess(:,:,5)=[0,0,-T2*sinTheta,0,0;
                              0,0,T2*cosTheta,0,0;
                              zeros(3,5)];
                if(nargout>3)
                    papt=zeros(5,1);
                end
            end
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
