function [xPred,MPred,DPred,PPred,Lambda]=reducedStateDiscPred(xPrev,MPrev,DPrev,F,G,LambdaChoice,modParams,u,du)
%%REDUCEDSTATEDISCPRED Perform the state prediction step of the discrete-
%               time or the discretized reduced state estimator. This is a
%               filter that takes measurements in Cartesian coordinates and
%               assumes that the dynamic model includes an input parameter
%               that depends on an unknown bounded parameter lambda. The
%               filter separates contributions due to measurement errors
%               and dynamic model mismatch errors. This implementation
%               allows for one to use general algorithms of [1] by
%               specifying F, G, a matrix for LambdaChoice, and Gammau,
%               Gammadu. However, by using the modParams input and omitting
%               other inputs, one can very easily use the maneuvering
%               models of [2] without having to explicitly create the
%               necessary inputs.
%
%INPUTS: xPrev The xDimX1 state estimate at the previous time-step.
%        MPrev The xDimXxDim matrix contributing to the total predicted
%              state covariance matrix based solely on measurement errors
%              at the previous time-step.
%        DPrev The xDimXzDim matrix of bias coefficients that are supposed
%              to relate target state errors to dynamic model parameter
%              uncertainty at the previous time-step.
%        F, G  F is the xDimXxDim state transition matrix for the dynamic
%              model xPred=F*xPrev+u, where u is an input. The matrix G is
%              an xDimXzDim matrix of partial derivatives of the control
%              input u with respect to the unknown bound parameter lambda.
%              If empty matrices are passed for these parameters, then xDim
%              must be a multiple of 2, and it is assumed that the state
%              consists of xDim/2 position components followed by xDim/2
%              velocity components and modParams.T (prediction time step)
%              is provided, so F=FPolyKal(2,xDim,1) and 
%              G=kron([T^2/2;T],eye(xDim/2,xDim/2)) which are the models
%              used for the filters in [2] for maneuvering models.
% LambdaChoice A parameter that specifies how the zDimXzDim (or 1X1)
%              covariance expressing the unknown parameter bounds is
%              formed. LambdaChoice can be the covariance matrix Lambda
%              itself (see Section IIB of [1] for one way this can be
%              formed). However, by setting the proper values in modParams
%              and passing a string for LambdaChoice, this function can
%              automatically compute Lambda for the two dynamic models in
%              [2]. Possible string values are
%                'GenTurn' A maximum linear acceleration A (typically
%                          meter/second^2) and turn rate Omega (radians/s)
%                          are provided in modParams. The filter forms
%                          Equation 8 of [2] using u based on the direction
%                          of the target in xPrev. Usually when this model
%                          is used, one will omit F and G to use the
%                          default dynamic values.
%               'GenAccel' A maximum acceleration is provided, but it is
%                          not specified in any particular direction
%                          (linear, tangential). This leads to the model
%                          for Lambda that is used in the simulations in
%                          [2]: Lambda=A*eye(xDim/2). Usually when this 
%                          model is used, one will omit F and G to use the
%                          default dynamic values.
%    modParams If F and G are not provied or if LambdaChoice is anything
%              but a matrix, then this is a structure containing the
%              parameters needed to use the default models. Possible
%              entries are
%              T The time that is used for the prediction step if F anf G
%                are not provied.
%              A The maximum linear acceleration if LambdaChoice='GenTurn'
%                and the maximum acceleration in any direction if
%                LambdaChoice='GenAccel'.
%              Omega The maximum turn rate (usually in radians per second)
%                if LambdaChoice='GenTurn'.
%       u, du u is an xDimX1 deterministic control input (which is a
%             function of the unknown parameter and for which the average
%             value of the unknown parameter is used) and du is a matrix of
%             its partial derivatives with respect to the elements of x. If
%             these parameters are omitted, default values of 0 are used.
%             Usually, when using LambdaChoice='GenTurn' or
%             LambdaChoice='GenAccel', these parameters will be omitted.
%
%OUTPUTS: xPred The xDimX1 predicted state estimate.
%         MPred The xDimXxDim predicted state covariance matrix
%               contribution due to measurement errors.
%         DPred The xDimXzDim matrix of predicted bias coefficients
%               contributing to the state covariance matrix.
%         PPred The xDimXxDim total covariance matrix of the predicted
%               state estimate.
%        Lambda The value of Lambda used in the state prediction.
%
%The algorithm is described generally in [1] and specifically to the
%default turning models used here in [2]. Unlike in [1], the control input
%here u is a full xDimX1 vector, not a reduced vector that is multiplied by
%some matrix to make it an xDimX1 vector.
%
%In [1] and [2], no clear method of initializing this type of tracking
%filter is provided. A simple way to initialize the filter would be to use
%two Cartesian converted measurements to obtain a state estimate and
%covariance as one would do with a normal Kalman filter (one could, for
%example, use the KalmanFIRSmoother function) and then set MPrev to the
%covariance value obtained while setting DPrev to zero.
%
%REFERENCES:
%[1] P. Mookerjee and F. Reifler, "Reduced state estimator for systems with
%    parametric inputs," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 2, pp. 446-461, Apr. 2004.
%[2] P. Mookerjee and F. Reifler, "Reduced state estimators for consistent
%    tracking of maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 41, no. 2, pp. 608-619, Apr. 2005.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
zDim=size(DPrev,2);

if(nargin<8||isempty(u))
    u=zeros(xDim,1);
end

if(nargin<9||isempty(Gammadu))
    du=zeros(xDim,xDim);
end

%If no dynamic model is specified, we assume the linear model used in [2]
%and assume that the user has provided T as an element of modParams. This
%only works if xDim is a multiple of 2.
if(isempty(F))
    T=modParams.T;
    F=FPolyKal(2,xDim,1);
    G=kron([T^2/2;T],eye(xDim/2,xDim/2));
end

if(~isa(LambdaChoice,'char'))
    %If the user provided a Lambda matrix.
    Lambda=LambdaChoice;
elseif(strcmp(LambdaChoice,'GenTurn'))%Maximum turn rate (and linear
                                      %acceleration) given.
    V=norm(xPrev((zDim+1):end));%The velocity
    
    %uVec is the unit vector in the direction of motion of the target (unit
    %tangent vector).
    uVec=xPrev((zDim+1):end)/V;
    if(any(~isfinite(uVec)))
        %This if-statement is to deal with the case where the target is not
        %moving:
        uVec=zeros(zDim,1);
        uVec(1)=1;
    end
    
    A=modParams.A;
    Omega=modParams.Omega;
    
    Lambda=A^2*(uVec*uVec')+V^2*Omega^2*(eye(zDim,zDim)-uVec*uVec');
elseif(strcmp(LambdaChoice,'GenAccel'))
    A=modParams.A;
    Lambda=A^2*eye(zDim,zDim);%Only a maximum acceleration is given.
else
    error('Unknown value given for lambdaChoice.');
end

FTotal=F+du;

%Equation 17 in [1].
MPred=FTotal*MPrev*FTotal';

%Ensure symmetry
MPred=(MPred+MPred')/2;

%Equation 18 in [1].
DPred=FTotal*DPrev+G;

%Equation 25 in [1].
xPred=F*xPrev+u;

%Equation 19 in [1].
PPred=MPred+DPred*Lambda*DPred';

%Ensure symmetry
PPred=(PPred+PPred')/2;
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
