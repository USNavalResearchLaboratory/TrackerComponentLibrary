function xSim=StrongStochRungeKStep(method,xCur,curT,a,D,deltaT,w,z,dadt,dDdt,dadx,d2adx)
%%STRONGSTOCHRUNGEKSTEP Simulate a nonlinear continuous-time random
%                       process specified by the Langevin equation forward
%                       in time by a step-size of deltaT using a strong
%                       Taylor scheme.
%
%INPUTS: method An integer defining which weak Taylor method to use.
%               1 The explicit order 1.5 strong Taylor scheme with additive
%                 noise. This requires the diffusion function to be
%                 additive.
%               2 The order 1.5 strong Taylor scheme with additive noise.
%                 This requires the diffusion function to be additive.
%          xCur The initial target state at time curT.
%          curT The time of the initial state xCur.
%             a The drift function in the continuous-time stochastic
%               dynamic model. It takes the state and a time variable as
%               its arguments, a(x,t).
%             D The diffusion function in the continuous-time stochastic
%               dynamic model. It takes the state and a time variable as
%               its arguments, D(x,t). If the process noise must be
%               additive, D(x,t) should not depend on x.
%        deltaT The size of the single step over which the strong
%               stochastic Runge-Kutta integration is performed.
%           w,z A pair of correlated normally distributed random variables
%               whose relationship is defined in section 10.4 of Kloeden.
%               If not provided or an empty matrix is passed, these values
%               are calculated. Either both or neither of these variables
%               must be supplied.
%          dadt A xDimX1 matrix of the derivative of the drift function
%               with respect to time at state xCur a   nd time curT. This
%               can also be a function that takes the state and a time
%               variable as its arguments, dadt(x,t).
%          dDdt A xDimXdColDim matrix of the derivative of the diffusion
%               function with respect to time at state xCur and time curT.
%               This can also be a function that takes the state and a
%               time variable as its arguments, dDdt(x,t). If the process
%               noise must be additive, dDdt(x,t) should not depend on x.
%          dadx A xDimXxDim matrix of the derivative of the drift function
%               with respect to the state at state xCur and time curT.
%               This can also be a function that takes the state and a
%               time variable as its arguments, dadx(x,t).
%         d2adx A xDimXxDimXxDim matrix of the second derivative of the
%               drift function with respect to the state at state xCur and
%               time curT. The value at point (m,k,l) represents
%               d2a(m)/dx(k)dx(l). This can also be a function that takes
%               the state and a time variable as its arguments,
%               d2adx(x,t). If not provided, this is assumed to be zero.
%
%OUTPUTS: xSim The simulated target state at time curT+deltaT.
%
%The algorithm for the explicit order 1.5 Taylor scheme is described in
%Eq. 2.19, Chapter 11.2 of [1]. The vector form of the equations used here
%is written out in [2]. The algorithm for the additive order 1.5 Taylor
%scheme is described in Chapter 10.4 of [1].
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, Part II, pp. 4-41, Feb. 2015.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

aCur=a(xCur,curT);
DCur=D(xCur,curT);
DNext=D(xCur,curT+deltaT);
dColDim=size(DCur,2);
xDim=size(aCur,1);

if(nargin<7 || isempty(w))
    %Generate the random component. It takes two correlated random variables.
    RSim=[deltaT*eye(dColDim),     deltaT^2/2*eye(dColDim);
        deltaT^2/2*eye(dColDim), deltaT^3/3*eye(dColDim)];
    SSim=chol(RSim,'lower');
    noiseVec=SSim*randn(dColDim*2,1);
    w=noiseVec(1:dColDim);
    z=noiseVec((dColDim+1):end);
end

switch(method)
    case 1 % Explicit order 1.5
        y=xCur+(1/dColDim)*aCur*deltaT;
        
        yp=bsxfun(@plus,y,DCur*sqrt(deltaT));
        ym=bsxfun(@minus,y,DCur*sqrt(deltaT));
        
        f1Val=xCur+deltaT*(1-dColDim/2)*aCur+(deltaT/4)*sum(a(yp,curT+deltaT)+a(ym,curT+deltaT),2);
        F3Val=1/(2*sqrt(deltaT))*(a(yp,curT+deltaT)-a(ym,curT+deltaT))-(1/deltaT)*(DNext-DCur);
        
        xSim=f1Val+DNext*w+F3Val*z;
    case 2 % Additive order 1.5
        if nargin<12
            d2adx=zeros(xDim,xDim,xDim);
        end
        
        if(isa(dadt,'function_handle'))
            dadt=dadt(xCur,curT);
        end
        if(isa(dDdt,'function_handle'))
            dDdt=dDdt(xCur,curT);
        end
        if(isa(dadx,'function_handle'))
            dadx=dadx(xCur,curT);
        end
        if(isa(d2adx,'function_handle'))
            d2adx=d2adx(xCur,curT);
        end
        
        %L Operators
        L0a=taylorOperator(aCur,DCur,0,dadx,dadt,d2adx);
        Lja=taylorOperator(aCur,DCur,[],dadx);
        
        %Additive order 1.5 strong Taylor scheme (Eq. 10.4.10 in Kloeden)
        f1Val=xCur+aCur*deltaT+(deltaT^2/2)*L0a;
        f2Val=DNext*w;
        f3Val=Lja*z+dDdt*(w*deltaT-z);
        
        xSim=f1Val+f2Val+f3Val;
    otherwise
        error('Not a valid method')
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
