function [aVal,aJacob,aHess,papt]=aCoordTurn2DTrans(xState,tauInv,tauAccelInv)
%%ACCORDTURN2DTRANS The continuous-time drift function for a 2D coordinated
%              turn model with a Cartesian state and a turn rate given as a
%              transversal acceleration. Additionally, a linear
%              acceleration can be given. The turn rate and linear
%              acceleration can optionally have time constants associated
%              with them, like in the Singer model, modelling a tendancy to
%              eventually want to return to non-accelerating, straight-line
%              motion. Note that the use of a transverse acceleration
%              component can be problematic when the velocity is near zero,
%              because the implied turn rate can diverge.  
%
%INPUTS: xState The 5X1 or 6X1 target state for 2D motion. If there is no
%              linear acceleration (acceleration along the direction of
%              motion), then xState=[x;y;xdot;ydot;at], where at is the
%              transversal acceleration. If there is a linear
%              acceleration, then the target state is
%              xState=[x;y;xdot;ydot;at;al], where al is the linear
%              acceleration. The dimensionality of the state is used to
%              determine whether a linear acceleration component is
%              present. The linear acceleration component changes the
%              speed. That means that it acts in the direction of the
%              velocity vector.
%       tauInv The inverse of the correlation time constant tau for the
%              transversal acceleration in seconds. tauInv must be
%              positive. The default if omitted or an empty matrix is
%              passed is zero.
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
%The concept of using the transversal acceleration instead of a turn rate
%is mentioned in [2], though no differential equations are given and a more
%detailed reference cited therein is a hard-to-get dissertation in French.
%The use of transversal acceleration is discussed in more detail in [3],
%though expressions are given when considering the 2D velocity are broken
%into components of heading and speed rather than in Cartesian space. The
%generalization to Cartesian space is not difficult and is done here.
%
%It can be shown that the magnitude of the transverse acceleration is equal
%to the speed times the turn rate. The model utilizing the transverse
%acceleration comes naturally from there. The relationship to the time
%constant in the transverse acceleration model comes from how it related to
%the derivative of the turn rate in [1].
%
%The optional linear acceleration provides a derivative of the speed of the
%target. The time constant operates in the same manner as the time constant
%for the turn rate. That is,
%alDot=-(1/tauLinAccel)*al+noise
%This is similar to how (total) acceleration decays in the Singer dynamic
%model, which was described in [4].
%
%The formulation in terms of transverse acceleration contains a singularity
%in the computation of the implied turn rate if the speed is zero. In such
%an instance, the turn rate is just set to zero. Similarly, singularities
%exist in applying the linear acceleration when the velocity of the target
%is zero. Values are zero are also substituted in such instances.
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
%[2] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[3] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design
%    of a multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%[4] R. A. Singer,"Estimating optimal tracking filter performance for
%    manned maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
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
    at=xState(5);%The transverse acceleration.
    xDot2=xDot*xDot;
    yDot2=yDot*yDot;

    v2=xDot2+yDot2;
    v=sqrt(v2);%The speed

    omega=at/v;
    if(~isfinite(omega))%Deal with zero velocity.
        zeroTerms=true;
        omega=0;
    else
        zeroTerms=false;
    end

    switch(length(xState))
        case 5%There is no linear acceleration.
        %This uses the relationship between turn rate, velocity and
        %transverse acceleration implied by Equation 2.A25 in
        %Blom's paper and puts it into Equations 61 and 67 in [1]. The
        %time-constant, which is not in [3], was added using the same logic
        %as [1], but substituting the appropriate terms in the expression
        %for the derivative of the turn rate.
            aVal=[xDot;%Position derivative
                  yDot;%Position derivative
                  -omega*yDot;%Velocity derivative
                  omega*xDot;%Velocity derivative
                  -tauInv*at];%Transverse acceleration derivative
            
            if(nargout>1)
                if(zeroTerms)
                    aJacob=[0,0,1,0,0;
                            0,0,0,1,0;
                            0,0,0,0,0;
                            0,0,0,0,0;
                            0,0,0,0,-tauInv];
                    if(nargout>2)
                        aHess=zeros(5,5,5); 
                        if(nargout>3)
                            papt=zeros(5,1);
                        end
                    end
                    return
                end

                v3=v*v2;

                pXpxDot=(at*xDot*yDot)/v3;
                pYpxDot=(at*yDot2)/v3;

                pXpyDot=-((at*xDot2)/v3);
                pYpyDot=-((at*xDot*yDot)/v3);

                pXpat=-yDot/v;
                pYpat=xDot/v;

                aJacob=[0,0,1,      0,      0;
                        0,0,0,      1,      0;
                        0,0,pXpxDot,pXpyDot,pXpat;
                        0,0,pYpxDot,pYpyDot,pYpat;
                        0,0,0,      0,      -tauInv];
                
                if(nargout>2)
                    aHess=zeros(5,5,5);
                    
                    v5=v3*v2;
                    
                    pXpxDotpxDot=(at*yDot*(-2*xDot2+yDot2))/v5;
                    pXpxDotpyDot=(at*xDot*(xDot2-2*yDot2))/v5;
                    pXpxDotpat=(xDot*yDot)/v3;
                    
                    pYpxDotpxDot=-((3*at*xDot*yDot2)/v5);
                    pYpxDotpyDot=-((at*yDot*(-2*xDot2+yDot2))/v5);
                    pYpxDotpat=yDot^2/v3;
                    
                    pXpyDotpxDot=pXpxDotpyDot;
                    pXpyDotpyDot=(3*at*xDot2*yDot)/v5;
                    pXpyDotpat=-(xDot^2/v3);
                    
                    pYpyDotpxDot=pYpxDotpyDot;
                    pYpyDotpyDot=-((at*xDot*(xDot2-2*yDot2))/v5);
                    pYpyDotpat=-((xDot*yDot)/v3);
                    
                    pXpatpxDot=pXpxDotpat;
                    pXpatpyDot=pXpyDotpat;
                    
                    pYpatpxDot=pYpxDotpat;
                    pYpatpyDot=pYpyDotpat;

                    aHess(:,:,3)=[0,0,0,            0,              0;
                                  0,0,0,            0,              0;
                                  0,0,pXpxDotpxDot, pXpyDotpxDot,   pXpatpxDot;
                                  0,0,pYpxDotpxDot, pYpyDotpxDot,   pYpatpxDot;
                                  0,0,0,            0,              0];

                    aHess(:,:,4)=[0,0,0,            0,              0;
                                  0,0,0,            0,              0;
                                  0,0,pXpxDotpyDot, pXpyDotpyDot,   pXpatpyDot;
                                  0,0,pYpxDotpyDot, pYpyDotpyDot,   pYpatpyDot;
                                  0,0,0,            0,              0];

                    aHess(:,:,5)=[0,0,0,            0,              0;
                                  0,0,0,            0,              0;
                                  0,0,pXpxDotpat,   pXpyDotpat,     0;
                                  0,0,pYpxDotpat,   pYpyDotpat,     0;
                                  0,0,0,            0,              0];
                    if(nargout>3)
                        papt=zeros(5,1);
                    end        
                end
            end
        case 6%There is a linear acceleration component.
            al=xState(6);
            %Get the linear acceleration vector. There will be NaNs if
            %the velocity is zero.
            linAccel=([xDot;yDot]/v)*al;

            %Deal with problems with a zero or nearly zero velocity.
            if(any(~isfinite(linAccel)))
                zeroTerms=true;
                linAccel(:)=0;
            end 

            aVal=[xDot;%Position derivative
                  yDot;%Position derivative
                  -omega*yDot+linAccel(1);%Velocity derivative
                  omega*xDot+linAccel(2);%Velocity derivative
                  -tauInv*at;%Transverse acceleration derivative
                  -tauAccelInv*al];%Linear acceleration derivative

            if(nargout>1)
                if(zeroTerms)
                    aJacob=[0,0,1,0,0,      0;
                            0,0,0,1,0,      0;
                            0,0,0,0,0,      0;
                            0,0,0,0,0,      0;
                            0,0,0,0,-tauInv,0;
                            0,0,0,0,0,      -tauAccelInv];

                    if(nargout>2)
                        aHess=zeros(6,6,6); 
                        if(nargout>3)
                            papt=zeros(6,1);
                        end
                    end

                    return;
                end

                v3=v*v2;
                pXpxDot=(yDot*(at*xDot+al*yDot))/v3;
                pXpyDot=-((xDot*(at*xDot+al*yDot))/v3);
                pXpat=-(yDot/v);
                pXpal=xDot/v;
                
                pYpxDot=(yDot*(-al*xDot+at*yDot))/v3;
                pYpyDot=(xDot*(al*xDot-at*yDot))/v3;
                pYpat=xDot/v;
                pYpal=yDot/v;
                
                aJacob=[0,0,1,      0,      0,      0;
                        0,0,0,      1,      0,      0;
                        0,0,pXpxDot,pXpyDot,pXpat,  pXpal;
                        0,0,pYpxDot,pYpyDot,pYpat,  pYpal;
                        0,0,0,      0,      -tauInv,0;
                        0,0,0,      0,      0,      -tauAccelInv];
                
                if(nargout>2)
                    aHess=zeros(6,6,6);
                    
                    v5=v3*v2;
                    xDot3=xDot2*xDot;
                    yDot3=yDot2*yDot;
                    
                    pXpxDotpxDot=(yDot*(-3*al*xDot*yDot+at*(-2*xDot2+yDot2)))/v5;
                    pXpxDotpyDot=(at*xDot3+2*al*xDot2*yDot-2*at*xDot*yDot2-al*yDot3)/v5;
                    pXpxDotpat=(xDot*yDot)/v3;
                    pXpxDotpal=yDot2/v3;
                    
                    pXpyDotpxDot=pXpxDotpyDot;
                    pXpyDotpyDot=(xDot*(-al*xDot2+3*at*xDot*yDot+2*al*yDot2))/v5;
                    pXpyDotpat=-(xDot2/v3);
                    pXpyDotpal=-((xDot*yDot)/v3);
                    
                    pXpatpxDot=pXpxDotpat;
                    pXpatpyDot=pXpyDotpat;
                    pXpatpat=0;
                    pXpatpal=0;
                    
                    pXpalpxDot=pXpxDotpal;
                    pXpalpyDot=pXpyDotpal;
                    pXpalpat=0;
                    pXpalpal=0;
                    
                    pYpxDotpxDot=-((yDot*(3*at*xDot*yDot+al*(-2*xDot2+yDot2)))/v5);
                    pYpxDotpyDot=(-al*xDot3+2*at*xDot2*yDot+2*al*xDot*yDot2-at*yDot3)/v5;
                    pYpxDotpat=yDot^2/v3;
                    pYpxDotpal=-((xDot*yDot)/v3);
                    
                    pYpyDotpxDot=pYpxDotpyDot;
                    pYpyDotpyDot=-((xDot*(3*al*xDot*yDot+at*(xDot2-2*yDot2)))/v5);
                    pYpyDotpat=-((xDot*yDot)/v3);
                    pYpyDotpal=xDot^2/v3;
                    
                    pYpatpxDot=pYpxDotpat;
                    pYpatpyDot=pYpyDotpat;
                    pYpatpat=0;
                    pYpatpal=0;
                    
                    pYpalpxDot=pYpxDotpal;
                    pYpalpyDot=pYpyDotpal;
                    pYpalpat=0;
                    pYpalpal=0;
                    
                    aHess(:,:,3)=[0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0;
                                  0,0,pXpxDotpxDot,pXpyDotpxDot,pXpatpxDot,  pXpalpxDot;
                                  0,0,pYpxDotpxDot,pYpyDotpxDot,pYpatpxDot,  pYpalpxDot;
                                  0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0];
                    aHess(:,:,4)=[0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0;
                                  0,0,pXpxDotpyDot,pXpyDotpyDot,pXpatpyDot,  pXpalpyDot;
                                  0,0,pYpxDotpyDot,pYpyDotpyDot,pYpatpyDot,  pYpalpyDot;
                                  0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0];
                    aHess(:,:,5)=[0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0;
                                  0,0,pXpxDotpat,pXpyDotpat,pXpatpat,  pXpalpat;
                                  0,0,pYpxDotpat,pYpyDotpat,pYpatpat,  pYpalpat;
                                  0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0];
                    aHess(:,:,6)=[0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0;
                                  0,0,pXpxDotpal,pXpyDotpal,pXpatpal,  pXpalpal;
                                  0,0,pYpxDotpal,pYpyDotpal,pYpatpal,  pYpalpal;
                                  0,0,0,      0,      0,      0;
                                  0,0,0,      0,      0,      0];
                              
                    if(nargout>3)
                        papt=zeros(5,1);
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
