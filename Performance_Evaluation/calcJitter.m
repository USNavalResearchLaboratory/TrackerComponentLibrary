function J=calcJitter(xEst,param2)
%%CALCJITTER Determine the "jitter" of a target track using one of two
%            definitions. One definitions deals with how the estimation
%            error varies when given truth data. The other defintion can be
%            used in the absence of truth data and deals with how well
%            predicted values agree with observed values.
%
%INPUTS: xEst This is a numDimXnumTargetsXnumSteps set of positions (if
%             truth is given) or target state esimates (if no truth is
%             given) whose jitter should be computed. Note that numSteps>1.
%      param2 The nature of this parameter determines the algorithm that is
%             used. This can be xTrue, a numDimXnumTargetsXnumSteps set of
%             true target position values, in which case the algorthm in
%             [1] using truth data is executed. Otherwise, this must be a
%             structure with members that predict the values in xEst
%             forward in time and extract the position components and the
%             algorithm of [1] that does not require truth data is used. To
%             predict forward in time, if the system is linear, then the
%             matrix param2.F must be provided. If the prediction is
%             nonlinear, then the function handle param2.f must be
%             provided. Similarly, for position components coming linearly
%             from the state, then the matrix param2.H must be provided.
%             Otherwise, the function handle param2.h must be provided.
%
%OUTPUTS: J A numTarX1 set of jitter measurements for each target.
%
%Note that when dealing with a single track, one must make sure that the
%matrix has a second dimension of 1.
%
%EXAMPLE:
%Here, we generate a random linear track as the truth. We then add noise to
%all of the estimates and compute the jitter using both of the methods.
% T=1;
% order=1;
% numDim=4;
% F=FPolyKal(T,numDim,order);
% q0=1;
% Q=QPolyKal(T,numDim,order,q0);
% SQ=chol(Q,'lower');
% %A matrix to extract position for a linear model.
% H=[1,0,0,0;
%    0,1,0,0];
% 
% numSteps=100;
% xInit=[0;0;100;0];
% numTargets=1;
% xTrue=zeros(numDim,numTargets,numSteps);
% xTrue(:,1,1)=xInit;
% for curStep=2:numSteps
%     xTrue(:,1,curStep)=F*xTrue(:,curStep-1)+SQ*randn(numDim,1);
% end
% 
% %We now have a true trajectory. The jitter with respect to itself is
% %zero (we only use position components).
% J0=calcJitter(xTrue(1:2,:,:),xTrue(1:2,:,:))
% %However, the jitter without the truth is not zero; it is driven by the
% %process noise.
% param2.H=H;
% param2.F=F;
% J1=calcJitter(xTrue,param2)
% %On the other hand, if we add noise both types of jitter will increase.
% 
% xNoisy=xTrue+10*rand(numDim,1,numSteps);
% J2=calcJitter(xNoisy(1:2,:,:),xTrue(1:2,:,:))
% J3=calcJitter(xNoisy,param2)
%The errors are both larger, though the error without truth is still
%notably larger than that with truth.
%
%REFERENCES:
%[1] O. E. Drummond, "Methodologies for performance evaluation of
%    multitarget multisensor tracking," in Proceedings of SPIE: Signal and
%    Data Processing of Small Targets Conference, vol. 3809, Jul. 1999, pp.
%    355-369.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTargets=size(xEst,2);
numSteps=size(xEst,3);

J=zeros(numTargets,1);

if(isnumeric(param2))
    xTrue=param2;
    xErr=xEst-xTrue;
    for curTarget=1:numTargets
        JCur=0;
        for curStep=2:numSteps
            diff=xErr(:,curTarget,curStep)-xErr(:,curTarget,curStep-1);
            JCur=JCur+diff'*diff;
        end
        J(curTarget)=JCur;
    end
else
    if(isfield(param2,'H')&&isfield(param2,'F'))
        H=param2.H;
        F=param2.F;

        for curTarget=1:numTargets
            JCur=0;
            for curStep=2:numSteps 
                diff=H*xEst(:,curTarget,curStep)-H*F*xEst(:,curTarget,curStep-1);
                JCur=JCur+diff'*diff;
            end
            J(curTarget)=JCur;
        end
    elseif(isfield(param2,'H'))
        H=param2.H;
        f=param2.f;
        for curTarget=1:numTargets
            JCur=0;
            for curStep=2:numSteps 
                diff=H*xEst(:,curTarget,curStep)-H*f(xEst(:,curTarget,curStep-1));
                JCur=JCur+diff'*diff;
            end
            J(curTarget)=JCur;
        end
    elseif(isfield(param2,'F'))
        h=param2.h;
        F=param2.F;

        for curTarget=1:numTargets
            JCur=0;
            for curStep=2:numSteps 
                diff=h(xEst(:,curTarget,curStep))-h(F*xEst(:,curTarget,curStep-1));
                JCur=JCur+diff'*diff;
            end
            J(curTarget)=JCur;
        end
    else
        h=param2.h;
        f=param2.f;

        for curTarget=1:numTargets
            JCur=0;
            for curStep=2:numSteps 
                diff=h(xEst(:,curTarget,curStep))-h(f(xEst(:,curTarget,curStep-1)));
                JCur=JCur+diff'*diff;
            end
            J(curTarget)=JCur;
        end
    end
end

J=sqrt(J);

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
