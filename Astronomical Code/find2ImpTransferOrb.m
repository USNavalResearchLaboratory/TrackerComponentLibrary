function [deltaV1,deltaV2,TOpt]=find2ImpTransferOrb(stateVec1,stateVec2,costSel,TMax,GM)
%%FIND2IMPTRANSFERORB Figure the optimal two-impulse maneuvers necessary to
%                 get a spacecraft from Cartesian position and velocity
%                 state vector stateVec1 to stateVec2. The impulse forces
%                 are applied at the beginning and at the end of the
%                 trajectory. Only a simpler Keplerian two-body system is
%                 considered. The maneuver is optimal either in that it
%                 minimizes the sum of the magnitudes of the change in the
%                 velocities (minimizing Newtonian kinetic energy/ fuel) or
%                 in minimizing the sum of the squares of the changes in
%                 the velocities. Only a zero-revolution solution is
%                 considered. Thus, this is not appropriate for Hohmann
%                 transfers as a much lower energy solution can be obtained
%                 by allowing for a single orbit to complete the transfer.
%
%INPUTS: stateVec1 A 6X1 vector consisting of position and velocity in a
%                  Cartesian inertial coordinate system where the
%                  gravitating body is at the origin.The units are assumed
%                  to be meters and meters per second.
%        stateVec2 The desired 6X1 state vector after the orbital maneuver.
%          costSel An optional string specifying the optimality criterion.
%                  Possible values are:
%                  'minAbsVelSum' (the default if omitted) Minimize the sum
%                                 of the magnitudes of the changes in
%                                 velocity for the impulse maneuvers.
%                  'minVelSquaredSum' Minimize the sum of the squares of
%                                 the mangnitudes of the changes in
%                                  velocity for the impulse maneuvers.
%             TMax An optional parameter setting an upper bound on the time
%                  of the optimal transfer. If omitted, the upper bound is
%                  set to 10 times the minimum time elliptical orbit
%                  between the points given by orbVelDet2PtMinEng. This
%                  parameter affects the convergence of the algorithm as a
%                  bounded line search method is used.
%               GM An optional value of the universal gravitational
%                  constant times the mass of the Earth. If omitted, the
%                  value Constants.WGS84GMWithAtmosphere is used. The
%                  units are km^3/sec^2.
%
%The function orbVelDet2Pt determines the velocities of an orbit given two
%points and a propagation time. The propagation time needed to optimize the
%selected cost function is unknown. This function solves for the optimal
%propagation time that minimizes the cost function using the fminbnd
%function, which performs a line search over a bounded region.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    r1Vec=stateVec1(1:3,1);
    V1Vec=stateVec1(4:6,1);
    r2Vec=stateVec2(1:3,1);
    V2Vec=stateVec2(4:6,1);
    
    if(nargin<5)
        GM=Constants.WGS84GMWithAtmosphere;
    end

    %TMin and TMax bound the search region.
    TMin=0;
    %A reasonable guess for an upper bound on the time that the minimum
    %energy transfer orbit would take is 10 times the maximum time that a
    %minimum energy orbit between the points would take.
    if(nargin<4)
       [~,TMax]=orbVelDet2PtMinEng(r1Vec,r2Vec,GM);
       TMax=10*TMax;
    end
    
    if(nargin<3)
        costSel='minAbsVelSum';
    end
    
    switch(costSel)
        case 'minAbsVelSum'
            costFun=@(T)costFunMinFuel(T);
        case 'minVelSquaredSum'
            costFun=@(T)costFunMinSquaredVel(T);
        otherwise
            error('An invalid cost function was specified.')
    end
    
    TOpt=fminbnd(costFun,TMin,TMax);
    [~,deltaV1,deltaV2]=costFun(TOpt);
    
    function [costVal,deltaV1,deltaV2]=costFunMinFuel(T)
        %Find the velocities of the transfer orbit.
        [W1,W2]=orbVelDet2Pt(r1Vec,r2Vec,T,0,false,false,GM);
        
        deltaV1=W1-V1Vec;
        deltaV2=W2-V2Vec;
        costVal=norm(deltaV1)+norm(deltaV2);
    end

    function [costVal,deltaV1,deltaV2]=costFunMinSquaredVel(T)
        %Find the velocities of the transfer orbit.
        [W1,W2]=orbVelDet2Pt(r1Vec,r2Vec,T,0,false,false,GM);
        
        deltaV1=W1-V1Vec;
        deltaV2=W2-V2Vec;
        costVal=norm(deltaV1)^2+norm(deltaV2)^2;
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
