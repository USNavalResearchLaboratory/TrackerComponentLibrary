function [knownSols,targetStates,exitFlag]=orbDet3Dir(numHalfRevs,obsLocs,unitDirVecs,t12,t13,solSel,findAll,knownSols,rho13EstVec,GM)
%%ORBDET3DIR Determine an orbital (ballistic) trajectory using 3 direction
%            vectors, the number of half-orbits of the satellite between
%            the first and third vectors, the time differences between the
%            first and second and between the first and third vectors,
%            and the location of the observer at each of the times. A
%            simple Keplerian dynamic model is assumes with the
%            mass-producing body at the origin. Multiple solutions to the
%            problem can exist and it is possible that no solutions might
%            exist. Additional solutions are obtained by calling the
%            function again with a new initialization and passing the
%            known solutions as a parameter, or by using the findAll
%            option for a particular set of numHalfRevs and solSel.
%            Not all problems have a solution. Invalid solutions 
%            (those with negative ranges from the observer) can be 
%            produced and must be passed as known solution back to the
%            function when searching for additional valid solutions. Not
%            all solutions will necessarily be found. Light-time is not
%            taken into account.
%
%INPUTS: numHalfRevs The number >=0 of half orbital revolutions between the
%                    first and third direction observation.
%            obsLocs A 3X3 matrix of the Cartesian location of the observer
%                    in a (quasi)-inertial coordinate system with the mass
%                    producing body at the origin at the times of each of
%                    the three direction observations. obsLocs(:,1) is the
%                    location at the time of the first observation
%                    (usually in meters). No observer should be at the
%                    origin.
%        unitDirVecs Unit direction vectors from the observer to the target
%                    at each of the three times in the global,
%                    (quasi)-inertial coordinate system. unitDirVecs(:,1)
%                    is the direction from the observer to the target at
%                    time 1 (unitless).
%                t12 The time interval between the first and second
%                    observations (>0), (usually in seconds).
%                t13 The time interval between the first and third
%                    observations (>0), (usually in seconds).
%             solSel If numHalfRevs>1, two possible sets of solutions can
%                    be found. This selects between the two solution sets.
%                    Possible values are 0 and 1. For numHalfRevs<=1, this
%                    should be zero or no solutions will be found. For a
%                    given problem, it is possible that there are no
%                    solutions in set 1 and solutions in set 2 or vice
%                    versa, so both sets should usually be checked. If this
%                    parameter is omitted, the default value of 0 is used.
%            findAll For a given set of numHalfRevs and solSel, find all of
%                    the solutions using the method specified by Gooding in
%                    his paper.
%          knownSols A 3XnumSols matrix of known solutions FOR A GIVEN
%                    VALUE OF numHalfRevs AND solSel, where knownSols(:,i)
%                    is the ith solution consisting of the range from the
%                    observer to the target at each of the three times in
%                    order. If there are no known solutions, then this
%                    parameter can be omitted or an empty matrix can be
%                    passed.
%        rho13EstVec A 2X1 vector holding initial estimates of the
%                    distances of the target at times 1 and 3. This affects
%                    to which solution the algorithm converges. If this
%                    parameter is omitted or an empty matrix is passed, a
%                    default initial guess based on the values of obsLocs
%                    is used.
%                 GM An optional value of the universal gravitational
%                    constant times the mass of the Earth. If omitted, the
%                    value Constants.WGS84GMWithAtmosphere is used. The
%                    units are usually m^3/sec^2.
%
%OUTPUTS: knownSols The input matrix knownSols augmented with the newly
%                   found solutions for the ranges from the observer to the
%                   target at each time. Invalid solutions (negative
%                   ranges) are possible and should be included when
%                   performing future calls of the function. A maximum of
%                   one solution will be found for each call to the
%                   function.
%      targetStates If a solution is found, this is a 6X3 matrix of the
%                   Cartesian position and velocity of the target at each
%                   of the three times under a Keplerian dynamic model. If
%                   findAll=true, then this is a 6X3XN set of the N
%                   solutions found. Otherwise, an empty matrix is
%                   returned if no solution is found.
%          exitFlag A number specifying the termination state of the
%                   algorithm. When no solutions exist, the algorithm will
%                   always end with a nonzero exit flag. Possible values
%                   are
%                  -1  The maximum number of iterations (100) was reached
%                      without convergence.
%                   0  A solution was found.
%                   1  The algorithm terminated because it could not find a
%                      valid estimate of the location of the target at time
%                      This is often because the Lambert problem (finding
%                      the target state at time 2 given the positions at
%                      times 1 and 3) could not be solved.
%                   2  A problem evaluating numerical dervatives arose.
%                      This can either be because no valid step size could
%                      be found, or because the gradient matrix used in
%                      Newton's method turned out singular and thus could
%                      not be inverted.
%
%The algorithm is based on the algorithm of Gooding, which is described in
%[1] and [2].
%
%Given initial estimates of the distance of the target from the observer at
%times 1 and times 3, the algorithm solves Lambert's problem using the
%function orbVelDet2Pt to determine the velocity of the target at time 1.
%The target's state (position and velocity) is then propagated to time 2
%using the function KeplerOrbitProp. If the observer-target distances at
%times 1 and 3 are correct, then the estimated target state at time 2 will
%be aligned with unitDirVecs(:,2). Gooding's algorithm correct for
%misalignment by using a combination of Halley's method and Newton's
%method. In this implementation, only Newton's method is used.
%
%To use Newtons' method, the cost function should have the same
%dimensionality as the state (the state being the distances to the target
%at time 1 and time 3) so that a unique solution for each step exists.
%Gooding uses cross products to set up an orthogonal coordinate system
%where unitDirVecs(:,2) and the current computed estimate of the target at
%time p2, p2Est, form a plane. The primary penalty function f is the
%projection of p2Est onto a vector in the plane orthogonal to
%unitDirVecs(:,2). If the two vectors align, then this projection will have
%zero magnitude. The second penalty function g keeps intermediate estimates
%in the plane. Thus, the second penalty function is the projection onto the
%normal vector of the plane. As p2Est defines the plane, this projection is
%zero. However, this definition is important for finding the derivatives
%necessary for Newton's method. The derivatives are found numerically using
%the numDiff function.
%
%The cost function to determine convergence is the primary penalty function
%f divided by max(cr,r2), where cr=dot(p2Est,unitDirVecs(:,2)) and r2 is
%the distance from the origin to the observer at time 2. The cost function
%is ad-hoc and follows from Gooding's suggestion in the text of his paper.
%
%The cost function is not convex and a step can worsen convergence.
%Following Gooding's suggestions, if a step more than doubles the cost,
%then the step size is reduced using a secant/ regula falsi method. The
%cost function is only reduced twice before another iterations is
%attempted. Forcing the cost function to always be reduced during a step
%appears to lessen the probability to finding additional solutions.
%
%Unlike in Gooding's work, an additional method of step size reduction is
%attempted if a step produces no solution. In such an instances, the step
%size is divided a few times before giving up and trying a new ad-hoc
%initialization. The ad-hoc initialization methods are from Gooding.
%Re-initialization is tried every time the algorithm fails until all ad-hoc
%initialization routines have been used.
%
%The algorithm is not a translation of the code given in the appendices of
%Gooding's paper and report, as the code was written in Fortran and uses a
%number of goto instructions, which are not supported in Matlab.
%
%When searching for multiple solutions,  the modified cost function
%described by Gooding is used to discourage convergence to old solutions.
%
%Additional information on the three-vector orbit determination problem is
%given in Chapter 7.3 of [3].
%
%The algorthm can be demonstrated using the data from the revised first
%example of Escobal in Gooding's report:
% GM=11467.89807;
% obsLocs=[0.16606957, -0.73815134, -0.73343987;
%          0.84119785, -0.41528280, -0.42352540;
%         -0.51291356,  0.53035336,  0.53037164];
% unitDirVecs=[-0.92475472, 0.80904274, 0.85044131;
%              -0.37382824,-0.55953385,-0.49106628;
%              -0.07128226, 0.17992142, 0.18868888];
% %Make sure that the unit vectors are truly normalized.
% unitDirVecs=bsxfun(@rdivide,unitDirVecs,sqrt(sum(unitDirVecs.*unitDirVecs,1)));
% t12=0.0381533;
% t13=0.0399364;
% numHalfRevs=0;
% solSel=0;
% %This will find all of the solutions that it can using the default
% %initialization/ trying others if that one fails.
% [knownSols,targetStates,exitFlag]=orbDet3Dir(numHalfRevs,obsLocs,unitDirVecs,t12,t13,solSel,true,[],[],GM)
%
%Four solutions should have been found. Specifically:
% knownSols =
%
%    1.100072659039377   3.591011260337803   0.014616956244466   7.508030332430836
%    1.518998485735830   1.883891729167911   0.446474185383019   2.970622569126413
%    1.578206974021349   2.031800780361934  -0.291261276964308   3.290782834714356
%
%Another example is the Mardsen example in Gooding's paper:
% GM=1;
% obsLocs=[0.976158959, -0.570816853, -0.909205804;
%         -0.228543158,  0.735996917,  0.359348664;
%         -0.099056153,  0.319119450,  0.155811658];
% unitDirVecs=[0.981950322,-0.515268519,-0.420980635;
%              0.176356395, 0.759734316, 0.802047494;
%             -0.068351936, 0.396613317, 0.423668647];
% t12=71.51143073;
% t13=72.04293220;
% %Make sure that the unit vectors are truly normalized.
% unitDirVecs=bsxfun(@rdivide,unitDirVecs,sqrt(sum(unitDirVecs.*unitDirVecs,1)));
% numHalfRevs=20;
% solSel=1;
% %This will find all of the solutions that it can using the default
% %initialization/ trying others if that one fails.
% [knownSols,targetStates,exitFlag]=orbDet3Dir(numHalfRevs,obsLocs,unitDirVecs,t12,t13,solSel,true,[],[],GM)
%
%This will find one solution, specifically,
% knownSols =
% 
%    0.054858143945831
%    0.066613740678771
%    0.095173924103870
%
%This is a 20 helf-revolution solution, so a solution coud potentiall have
%existed with solSel=0. However, had that been tried in this example, no
%solutions would have been found.
%
%REFERENCES:
%[1] R. H. Gooding, "A new procedure for the solution of the classical
%    problem of minimal orbit determination from three lines of sight,"
%    Celestial Mechanics and Dynamical Astronomy, vol. 66, no. 4, pp.
%    387-423, 1997.
%[2] R. H. Gooding, "A new procedure for orbit determination based on three
%    lines of sight (angles only)," Royal Aerospace Executive, Procurement
%    Executive, Ministry of Defence, Farnborough, Hants, United Kingdom,
%    Tech. Rep. 93004, Apr. 1993.
%[3] D. A. Vallado and W. D. McClain, Fundamentals of Astrodynamics and
%    Applications, 4th ed. Hawthorne, CA: Microcosm press, 2013.
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10)
    GM=Constants.WGS84GMWithAtmosphere;
end
  
%If no initial estimate is given, start with the ad-hoc initialization
%suggested by Gooding.
if(nargin<9||isempty(rho13EstVec))
    initMag=2*sum(norm(obsLocs(:,3))+norm(obsLocs(:,2))+norm(obsLocs(:,1)));
    rho13EstVec=[initMag;initMag];
end
 
%If no known solutions are provided.
if(nargin<8)
    knownSols=[];
    numKnownSol=0;
    r13SqKnownObs=[];
elseif(~isempty(knownSols))
    numKnownSol=size(knownSols,2);
    %This is to hold the squared distances to the known solutions of the 
    %target location at times t1 and t3 for use in the penalized cost
    %function that is minimized.
    r13SqKnownObs=zeros(2,numKnownSol);
    for curKnownSol=1:numKnownSol
        r13SqKnownObs(1,curKnownSol)=(norm(obsLocs(:,1)+knownSols(1,curKnownSol)*unitDirVecs(:,1)))^2;
        r13SqKnownObs(2,curKnownSol)=(norm(obsLocs(:,3)+knownSols(3,curKnownSol)*unitDirVecs(:,3)))^2;
    end
else
    r13SqKnownObs=[];
    numKnownSol=0;
end
 
if(nargin<7)
    findAll=false;
end
 
%If the solution branch is not specified.
if(nargin<6)
    solSel=0;
end
 
if(findAll)
    %If findAll is set, this will start with the given initialization, if
    %any, and then try to find all of the solutions for the given
    %combination of solSel and numHalfRevs using the ad-hoc
    %reinitialization method of Gooding.
    findAll=false;
    
    %Get an initial solution
    [knownSols,targetStates,exitFlag]=orbDet3Dir(numHalfRevs,obsLocs,unitDirVecs,t12,t13,solSel,findAll,knownSols,rho13EstVec,GM);
    if(exitFlag==0)%If a solution was found
        %Get a new estimate as Gooding described in Section 3.5 of his
        %report.
        rho13EstVec=2*knownSols([1,3])-rho13EstVec;
        curSol=2;
        while(1)
            %Try to get a new set of solutions.
            [knownSols,curTargetStates,exitFlag]=orbDet3Dir(numHalfRevs,obsLocs,unitDirVecs,t12,t13,solSel,findAll,knownSols,rho13EstVec,GM);
            if(exitFlag~=0)
                break;
            end
            %Augment the target states with the newly solved state
            %solution.
            targetStates(:,:,curSol)=curTargetStates;
            rho13EstVec=2*knownSols([1,3],1)-knownSols([1,3],2);
            curSol=curSol+1;
        end
    end
    return 
end

%The precision for declaring convergence.
convergeVal=1e-12;
maxIter=100;
 
%This is the fraction of the estimated origin-target range used for a
%step for numerical differentiation.
numDiffStepFrac=1e-5;
%The maximum number of times the step size for numerical differentiation is
%reduced before declaring the problem a failure.
maxDiffStepReductions=5;
 
%If a step does not reduce the cost function enough, then the secant/
%regula falsi method described by Gooding is used. maxNumBisections is the
%number of times such a reduction is performed.
maxNumBisections=2;
%On the other hand, if no solution to Lambert's function can be found for a
%given step size, then this is the number of attempts that will be made at
%reducing the step size before giving up.
maxNumStepReductions=5;
 
%This is the number of times ad-hoc reintialization routines have been
%attempted. Once all ad-hoc reinitialization routines fail, the algorithm
%is terminated.
curInitNumber=0;
%The number of reinitialization routines coded in this function. It will
%try all initialization approahces before declaring a total failure.
maxInits=3;
 
%Setting these to empty matrices tells the getP2ForNewIter function that
%no previous solution exists.
rho13OldVec=[];
fCostOld=[];
algConverged=false;
 
%Get the estimated location of the target at time t2 with respect to
%the observer p2Est, also get fc, the cost function for determining
%convergence and fCost, the cost function over which optimization is
%performed. rho13EstVec is returned again along with curInitNumber in the
%event that a failure occurred and a new intiialization was tried. If p2Est
%is empty, then a total failure occurred. Only parameters that are modified
%are passed to the function; the function implicitly takes a number of
%other parameters.
[p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=getP2ForNewIter(rho13EstVec,curInitNumber);
 
%Terminate with failure if no valid initial estimate could be found.
if(isempty(p2Est))
    exitFlag=1;
    return;
end

curIter=1;
while(1)    
    %Perform first-order numerical differentiation of the cost function.
    %If it fails, it is done a few more times with smaller step sizes for
    %the numerical differentiation. The initial step size is set to
    %[numDiffStepFrac times the distance of target 1 (for rho1Est) and of
    %target 3 from the global origin.
    p1=norm(obsLocs(:,1)+rho13EstVec(1)*unitDirVecs(:,1));
    p3=norm(obsLocs(:,3)+rho13EstVec(2)*unitDirVecs(:,3));
 
    %Loop until a good numerical derivative is obtained or the maximum
    %number of reductions of the differential step size is exceeded.
    numDiffStepReductions=0;
    prevInitNumber=curInitNumber;
    while(1)
        stepSize=max(numDiffStepFrac*[p1;p3]);
        %Second-order numerical differentiation.
        costFunDerivs=numDiff(rho13EstVec,@costFun,2,1,stepSize);        
        
        %If the derivatives were successfully computed.
        if(~isempty(costFunDerivs)&&~any(any(~isfinite(costFunDerivs))))
            break;
        end
        
        numDiffStepReductions=numDiffStepReductions+1;
        %Iff the number of step reductions exceeds the maximum number of
        %reductions, try reinitializing. If that fails, then exit.
        if(numDiffStepReductions>maxDiffStepReductions)
            %Get a new initialization.
            [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=reInit(curInitNumber);
            %If no new initialization could be found.
            if(isempty(p2Est))
                targetStates=[];
                exitFlag=2;
                return;
            else%Restart iterations on the new initalization
                curIter=1;
                break
            end
        end
        
        %If the numerical differentiation failed, try reducing the step
        %size.
        numDiffStepFrac=numDiffStepFrac/10;
    end
    %If the derivative routine reset the initialization, then start the
    %loop over.
    if(prevInitNumber~=curInitNumber)
        continue;
    end
    
    %Newton's method has to invert the derivative matrix. Thus, if it is
    %singular, the algorithm will not work. We check for that and indicate
    %a failure if the derivative matrix is singular. The extra check for
    %convergence of the cost function takes care of the case where the
    %derivative matrix is singular because the function converged --in
    %reality, that should only occur if the initial estimate is very close
    %to the truth.
    cr=dot(p2Est,unitDirVecs(:,2));
    r2=norm(obsLocs(:,2));
    if(rank(costFunDerivs)<2&&fc>convergeVal)
        %Get a new initialization.
        [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=reInit(curInitNumber);
        %If no new initialization could be found.
        if(isempty(p2Est))
            targetStates=[];
            exitFlag=2;
            return;
        else%Restart iterations on the new initalization
            curIter=1;
            continue;
        end
    elseif(rank(costFunDerivs)<2&&abs(fc/max(cr,r2))<=convergeVal)
        algConverged=true;
        break;
    end
 
    %When here, we have the numerical derivatives and can perform a step of
    %Newton's algorithm. The matrix costFunDerivs has the format:
    %costFunDerivs=[dfdrho1, dfdrho2
    %               dgdrho1, dgdrho2];
    %Since the cost function is not necessarily convex, there is no 100%
    %guarantee that the step will decrease the cost. Thus, after performing
    %the step, we must check that it has not worsened the cost function by
    %too much. The check and adjustment is done in the getP2ForNewIter
    %function.
 
    %Equation 6 in the report to get the Newton step.
    deltaRhoVec=-costFunDerivs\[fCost;0];
 
    fCostOld=fCost;
    rho13OldVec=rho13EstVec;
    rho13EstVec=rho13EstVec+deltaRhoVec;
    prevInitNumber=curInitNumber;
    [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=getP2ForNewIter(rho13EstVec,curInitNumber);
    %It it reinitialized, reset the iteration count.
    if(prevInitNumber~=curInitNumber)
        curIter=0;
    end
    
    %Terminate with failure if no current estimate of the location of the
    %target at time t2 could be found or at least none could be found that
    %reduces the step size.
    if(isempty(p2Est))
        targetStates=[];
        exitFlag=1;
        return;
    end
    
    %Prior to performing another iteration, we will check to see whether it
    %has already converged. The goal is to reduce fc to zero across
    %iterations. The convergence criterion is not entirely scale-invariant.
    cr=dot(p2Est,unitDirVecs(:,2));
    r2=norm(obsLocs(:,2));
    if(abs(fc/max(cr,r2))<convergeVal)
        algConverged=true;
        break;
    end
    curIter=curIter+1;
    
    if(curIter>maxIter)
        %If the maximum number of iterations has been exceeded without
        %convergence, then try reinitializing.
        [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=reInit(curInitNumber);
        %If no new initialization could be found.
        if(isempty(p2Est))
            break;
        else%Otherwise, start iterating the new initialization.
            curIter=1;
            continue;
        end
    end
end
 
%If we got here, then either the maximum number of iterations was exceeded
%or the algorithm converged.
if(algConverged)
    knownSols=[[rho13EstVec(1);norm(p2Est);rho13EstVec(2)],knownSols];
    exitFlag=0;
    return
else
    targetStates=[];
    exitFlag=-1;
    return
end
 
%The two-dimensional cost function
function costVec=costFun(rhoVec)
    p2Cur=findP2(numHalfRevs,obsLocs,unitDirVecs,rhoVec,t12,t13,GM,solSel);
    
    %If the findP2 function failed, return an empty matrix to indicate
    %failure of the derivative.
    if(isempty(p2Cur))
        costVec=[];
        return
    end
    
    %f and g are just projections onto the orthonormal vectors set up by
    %the last call to the getP2ForNewIter function.
    f=dot(eInPlane,p2Cur);
    g=dot(eNormal,p2Cur);
    
    %Adjust the cost functions so that known solutions are not repeated.
    if(~isempty(knownSols))
        for curKSol=1:numKnownSol
            xi=rhoVec(1)-knownSols(1,curKSol);
            eta=rhoVec(2)-knownSols(3,curKSol);
            beta=sqrt(xi^2+eta^2);
            gamma=sqrt(beta^2+sum(r13SqKnownObs(:,curKSol)));
            
            f=f*(gamma/beta);
            g=g*(gamma/beta);
        end
    end
    
    costVec=[f;g];
end
 
function [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=getP2ForNewIter(rho13EstVec,curInitNumber)
%To begin an iteration, the estimated location of the target at time t2
%with respect to the location of the observer at time t2 (obsLocs(:,2))
%implied by the current estimates of the distances from the observer to
%the target at times t1 and t2, and rho13EstVec, is needed. This
%current estimate will define the coordinate system used for the cost
%function when using Newton's method to update the step. However, if
%rho13EstVec is rather bad, then it is possible that no solutions exist.
%Thus, if the first attempt to find a solution fails, a few ad-hoc
%adjustments to rho1Est and rho3Est are undertaken. If these fail, then the
%algorithm is terminated as being unsuccessful.
%
%On the other hand, it it also possible that the algorithm might have been
%running for a few iterations with some initial estimate and then gotten
%stuck. In such an instance, the problem is reset and is treated as if the
%function had just begin. However, one does notwant to reset to the same
%value after it fails. Thus, the curInitNumber lists the last ad-hoc
%initialization tried before the function found a previous solution. Thus,
%the first time this function is run with whatever values the user
%supplied, curInitNumber=0 is a reasonable choice. If this function fails
%completely, then p2Est is returned empty.
%
%The function also returns fc, the cost value for convergence, and fCost,
%the cost value over which minization is performed.  Also, eNormal and
%eInPlane, which are necessary for evaluating the cost function, are
%returned.
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 
    [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
    
    %If an estimate is available, then we can go on to iteratively
    %improve out estimates of rho1Est and rho3Est as long as the estimate
    %is not worse than on the previous step. That is, if fCost<fCostOld.
    %Otherwise, we will try reducing the step size. The check for rho1Old
    %being empty takes care of the possibility that this is the first
    %function call.
    if(~isempty(p2Est)&&(isempty(rho13OldVec)||(abs(fCost)<2*abs(fCostOld))))
        return;
    end
    %If the attempt at getting a current solution failed, then some
    %attempt should be made to resurrect it. How this is done will
    %depend on whether or not this is the first iteration. If this is
    %not the first iteration, then no values for rho1Old, rho3Old,
    %or deltaRhoVec will be available, so we should just skip directly to
    %the ad-hoc initialization. Otherwise, we will first try reducing the
    %step size from the last iteration a few times.
    
    %If this is not the first iteration, then try reducing the step size
    %from the previous iteration. 
    if(~isempty(rho13OldVec))
        %This is the seacnt step size reduction suggested by Gooding
        if(~isempty(fCost))
            for curRed=1:maxNumBisections
                rho13EstVec=(fCostOld*rho13EstVec+rho13OldVec*fCost)/(fCostOld+fCost);
                
                [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
 
                if(~isempty(p2Est)&&abs(fCost)<2*abs(fCostOld))
                    return;
                elseif(isempty(fCost))
                    break;
                end
            end
        end
        if(isempty(fCost))
            %This step-size reduction method is used when the step pushes
            %the estimate into a region with no solution as the step size
            %reduction method of Gooding cannot be used without a cost at
            %the step location.
            for curRed=1:maxNumStepReductions
                rho13EstVec=rho13OldVec+deltaRhoVec/3^curRed;
                [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
 
                if(~isempty(p2Est))
                    return;
                end
            end
        else
            return;
        end
 
        %If no solution was found by reducing the step size, then we will
        %essentially reset the iterative process using some ad-hoc values.
    end
        
    %Try ad-hoc initialization routines.
    [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=reInit(curInitNumber);
        
    %If we got here, then the ad-hoc initialization routines failed. Return
    %an empty matrix for the initialization number to indicate failure.
end

function [p2Est,rho13EstVec,fc,fCost,curInitNumber,eNormal,eInPlane,targetStates]=reInit(curInitNumber)
%Try to reinitialize using an ad-hoc method.

    %Try ad-hoc initialization routines.
    while(curInitNumber<maxInits)
        if(curInitNumber==0)
        %The first ad-hoc initialization attempt is to use the projections
        %of the locations of the observers onto their line of sight lines,
        %as was done by Gooding.
            curInitNumber=curInitNumber+1;
            rho13EstVec(1)=-dot(unitDirVecs(:,1),obsLocs(:,1));
            rho13EstVec(2)=-dot(unitDirVecs(:,3),obsLocs(:,3));
            [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
            if(~isempty(p2Est))
                return;
            end
            %If it failed, we will try again with another ad-hoc method.
        end
        if(curInitNumber==1)
        %The second ad-hoc method is also by Gooding and involves
        %a series of projections of the locations of the observer with
        %respect to line-of-sight-values. He says that the result is the 
        %common perpendicular (shortest traversal) to the two lines of
        %sight.
            curInitNumber=curInitNumber+1;
            deltaObsLocs=obsLocs(:,3)-obsLocs(:,1);
            d1=dot(deltaObsLocs,unitDirVecs(:,1));
            d2=dot(unitDirVecs(:,1),unitDirVecs(:,3));
            d3=dot(deltaObsLocs,unitDirVecs(:,3));
            d4=1-d2^2;
            rho13EstVec(1)=(d1-d3*d2)/d4;
            rho13EstVec(2)=(d1*d2-d3)/d4;
            
            if(all(isfinite(rho13EstVec)))
                [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
                if(~isempty(p2Est))
                    return;
                end
            end
            %If it failed, we will try again with another ad-hoc method.
        end
        %The final ad-hoc technique is to just set the distance values to
        %zero.
        curInitNumber=curInitNumber+1;
        rho13EstVec=[0;0];
        [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols);
        if(~isempty(p2Est))
            return;
        end
    end
    
    %If here, initialization failed; just return empty matrices.
    p2Est=[];
    rho13EstVec=[];
    fc=[];
    fCost=[];
    curInitNumber=[];
    eNormal=[];
    eInPlane=[];
    targetStates=[];
end

end


 
 
function [p2Est,fc,fCost,eNormal,eInPlane,targetStates]=findP2(numHalfRevs,obsLocs,unitDirVecs,rho13EstVec,t12,t13,GM,solSel,r13SqKnownObs,knownSols)
%FINDP2 Estimate the Cartesian location of location of the satellite as
%       seen by the observer at obsLocs(:,2) at time t2 given estimates of
%       the ranges of the satellite at times t1 and t3. The estimate is
%       obtained by solving Lambert's problem using the location of the
%       satellite at times t1 and t3 (based on the direction angles and
%       range estimates) and propagating the estimate from time t1 forward
%       to time t2. The resulting location in global coordinates is then
%       translated to the coordinate system of the receiver at time t2.
%
%In the multirevolution case, it is possible for there to be two solutions.
%solSel selects which of the two solutions is to be used. solSel=0 chooses
%the first solution. solSel=1 chooses the second solution. However, if only
%one solution was found when solSel=1, then it returns as a failure (p2Est
%is empty).
%
%The function also returns the  the cost for convergence, fc, the cost over
%which iteration is performed, fCost, as well as eNormal and eInPlane,
%which are necessary for defining the cost function for numerical
%differentiation this iteration, and the states of the targets
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 
%These empty values will only be returned if the algorithm fails.
p2Est=[];
fc=[];
fCost=[];
targetStates=[];
eNormal=[];
eInPlane=[];
 
%First, find the Cartesian locations of the target at times 1 and 3 in
%global (ECI) coordinates.
tarLoc1=obsLocs(:,1)+rho13EstVec(1)*unitDirVecs(:,1);
tarLoc3=obsLocs(:,3)+rho13EstVec(2)*unitDirVecs(:,3);
 
%If either of the estimated points is at the global origin, then the
%algorithm will not work and we must return marking it as having failed.
if(all(tarLoc1==0)||all(tarLoc3==0))
    return;
end
 
%Find the total angle between the vctors, including the possibility of
%multiple revolutions between the points. For this, the numHalfRevs
%parameter is needed. It also helps disambiguate whether the long way or
%the short way in the orbit between the points is taken.
deltaTheta=angBetweenVecs(tarLoc1,tarLoc3);
%This switches to the long way around the orbit if an odd number of
%half-revolutions is given
if(mod(numHalfRevs,2))
    deltaTheta=pi-deltaTheta;
end
deltaTheta=deltaTheta+numHalfRevs*pi;
 
%Next, solve Lambert's problem to get the velocity of the target at point
%1.
[tarVel1,tarVel3]=orbVelDet2Pt(norm(tarLoc1),norm(tarLoc3),t13,deltaTheta,false,[],GM);
%There can be 0, 1, or two solutions.
numSol=size(tarVel1,2);
 
%If there are no solutions or if there is one solution, but solSel is two,
%indicating that we should use the second (nonexistent) solution, then
%return as a failure (p2Est is empty).
if(numSol==0||(numSol==1&&solSel==1))
    return;
elseif(solSel==1)%Use the second solution
    tarVel1=tarVel1(:,2);
    tarVel3=tarVel3(:,2);
else%Use the first solution
    tarVel1=tarVel1(:,1);
    tarVel3=tarVel3(:,1);
end
 
%The velocity tarVel1 is given parameterized in terms of a component in the
%direction of tarLoc1 (radial) and one in  the plane formed by tarLoc1
%and tarLoc3. We will convert it into global Cartesian coordinates.
uR1=tarLoc1/norm(tarLoc1);%Unit radial vector in the direction of tarLoc1.
 
%A vector normal to the plane in which tarLoc1 and tarLoc3 reside.
rNormal=cross(tarLoc1,tarLoc3);
%Normalize the vector.
uNormal=rNormal/norm(rNormal);
%A Unit vector in the tangent plane
uT1=cross(uNormal,uR1);
 
%If any of the values in the tangentional unit vectors is not finite, then
%the vectors are parallel (the trajectory is radial) and the tangential
%components of tarVel1 should be zero anyway, so just put in
%the proper values.
if(any(~isfinite(uT1)))
    tarVel1Cart=tarVel1(1)*uR1;
else%If the trajectory is not radial.
    tarVel1Cart=tarVel1(1)*uR1+tarVel1(2)*uT1;
end
targetState1=[tarLoc1;tarVel1Cart];
 
 
%The target state at point 1 is [tarLoc1;tarVel1Cart]. This is propagated
%by time interval t12 to get the estimated target state at t2 in global
%Cartesian coordinates
tarState2=KeplerOrbitProp(targetState1,t12,GM);
tarLoc2=tarState2(1:3,1);%The location in global Cartesian coordinates.
 
%The estimates location of the target at time t2 with respect to the
%observer's location at time t2.
p2Est=tarLoc2-obsLocs(:,2);
 
%If more than just p2Est for computing a cost is desired.
if(nargout>1)
    %p2Est and unitDirVecs(:,2) form a plane but are generally NOT orthogonal.
    %We want an orthogonal basis in which p2Est and unitDirVecs(:,2) are
    %parameterized. First, we get a vector perpendicular to the plane, and thus
    %perpendicular to unitDirVecs(:,2).
    eNormal=cross(unitDirVecs(:,2),p2Est);
    eNormal=eNormal/norm(eNormal);%Normalize
    %Now, we get a vector In the plane that is perpendicular to the normal and
    %to unitDirVecs(:,2) so that we have an orthogonal basis in the plane
    %itself.
    eInPlane=cross(eNormal,unitDirVecs(:,2));
    eInPlaneMag=norm(eInPlane);
    %Gooding performs the normalization every time f is computed; here, we do
    %it directly to the in-plane vector.
    eInPlane=eInPlane/eInPlaneMag;
 
    %As p2Est aligns with unitDirVecs(:,2), this goes to zero.
    fc=dot(eInPlane,p2Est);
    %Note that g=dot(eNormal,p2Est); is zero (within finite precision
    %limitations).
 
    fCost=fc;
    %Adjust the cost functions so that known solutions are not repeated. This
    %is taken from Section 3.5 of Gooding's report. However, instead of having
    %to save the full distances to the target, the sum of the distance
    %from the origin to the observer and the observer to the known target
    %location is used, which will always be greater than or equal to the
    %desired value.
    if(~isempty(knownSols))
        numKnownSols=size(knownSols,2);
        for curSol=1:numKnownSols
            xi=rho13EstVec(1)-knownSols(1,curSol);
            eta=rho13EstVec(2)-knownSols(3,curSol);
            beta=sqrt(xi^2+eta^2);
 
            %Note that one cannot just use the knownSols values, one must add
            %in something else or else if one has zero solutions, then
            %gamma=beta and the numerator and denominator cancel.
            gamma=sqrt(beta^2+sum(r13SqKnownObs(:,curSol)));
            fCost=fCost*(gamma/beta);
        end
    end
end
 
%The following is executed only if the matrix of target states is
%requested. Thus, it is not executed in the cost fucntion.
if(nargout>5)
    uR3=tarLoc3/norm(tarLoc3);
    uT3=cross(uNormal,uR3);
    
    if(any(~isfinite(uT1)))
        tarVel3Cart=tarVel3(1)*uR3;
    else%If the trajectory is not radial.
        tarVel3Cart=tarVel3(1)*uR3+tarVel3(2)*uT3;
    end
    
    targetState3=[tarLoc3;tarVel3Cart];
    targetStates=[targetState1,tarState2,targetState3];
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
