function [xCur,didConverge]=implicitEulerStep(x,curT,f,deltaT,givenVals,useNewton,maxIter,theta,RelTol,AbsTol)
%%IMPLICITEULERSTEP Perform a step of the implicit Euler method (also known
%                   as the backward Euler method) to integrate an ordinary
%                   differential equation (ODE) one step forward. Unlike
%                   the forward Euler method, the backward Euler method is
%                   L-stable. The function also takes a parameter alphaVal
%                   that can be varied to allow for others in the family of
%                   implicit Euler algorithms, such as the trapezoidal
%                   method. A fixed number of iterations is the default,
%                   but the algorithm can also be run until a convergence
%                   criterion is met.
%
%INPUTS: x The value of the xDimX1 state over which integration is being
%          performed.
%     curT The time t at which x is taken.
%        f A function handle such that [fCur,df]=f(x,t) returns the
%          derivative of x with respect to time taken at time t (fCur). If
%          one wishes to use Newton's method for the implicit iteration
%          (set useNewton=true), rather than just fixed point iteration,
%          then this function has to return a second value, being the
%          Jacobian of f, where each column is the derivative vector of f
%          with respect to the corresponding element of x. It is assumed
%          that the Jacobian is invertible; it is not checked for positive
%          definiteness.
%   deltaT The positive scalar size of the single step over which the
%          implicit step is taken. 
% givenVals This is a structure that provides fCur and df evaluated at x,
%          if available. The members of the structure are named fCur and
%          df. If useNewton=false, then fCur will be used. If
%          isNewton=true, then fCur and df will be used and if df is not
%          provided, then any value in givenVals is ignored and f is called
%          for both values. This input helps avoid repeatedly calculating
%          values of f that have already been found.
% useNewton If false, fixed point iteration is performed to take the step.
%          Otherwise, Newton's method is used (without a line search) and f
%          must return a Jacobian matrix as a second output. If omitted or
%          an empty matrix is passed, the default of false is used.
%  maxIter The maximum number of iterations to perform. If RelTol=0 and
%          AbsTol=0, then this equals the number of iterations performed.
%          The default if omitted or an empty matrix is 2. Setting this to
%          0 just results in the explicit Euler's method being used.
%    theta The parameter adjusting the effect of the prior on the
%          integration (see below). The default if omitted or an empty
%          matrix is passed is 1 (just the implicit Euler method. The
%          parameter varies from 0 to 1. A value of (1/2) makes this the
%          trapezoidal method.
% RelTol, AbsTol The relative and absolute tolerances on the iterations
%          before declaring convergence. If these are set to 0 (the default
%          if omitted or empty matrices are passed), then the algorithm
%          will just iterate for the maximum number of iterations. the
%          tolerances apply to each element of x. Convergence is declared
%          if all(diff<=AbsTol)||all(diff<=RelTol*abs(xCur)). Generally,
%          one will
%
%OUTPUTS: xCur The xDimX1 value of the x after taking an implicit Euler
%              step.
%  didConverge If RelTol and/or AbsTol are not zero and maxIter>0, then
%              this indicates whether the iterations converged to the
%              desired accuracy. Otherwise, this is just an empty matrix.
%
%The implicit Euler method, also known as the backward Euler method, which
%is also a first-order Adams-Moulton formula, solves the relation:
% x(curT+deltaT)=x(curT)+deltaT*f(x(curT+deltaT),curT+deltaT)
%Note that this differs from the forward Euler method, which is just.
% x(curT+deltaT)=x(curT)+deltaT*f(x(curT),curT)
%The forward Euler method comes from a truncated Taylor series expansion at
%curT and the backward Euler method can be derived from a truncated Taylor
%series expansion about curT+deltaT, where the point beign evaluated is
%curT.
%
%The implicit Euler method performs better with stiff equations than the
%forward Euler method, with the Newton iteration often the best choice for
%stiffer equations, since it converges faster.
%
%The family of implcit Euler methods allowed by the theta term have to do
%with replacing the above relation with 
%x(curT+deltaT)=x(curT)+deltaT*(theta*f(x(curT+deltaT),curT+deltaT)+(1-theta)*f(x(curT),curT))
%When theta=0.5, this is the trapezoidal rule, which is mentioned in
%Section 225 of Chapter 22 of [1].
%
%The implicit Euler method is discussed in Chapters 20 and 21 of [1].
%
%REFERENCES:
%[1] J. C. Butcher, Numerical Methods for Ordinary Differential Equations,
%    John Wiley and Sons: Chichester, 2003.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10||isempty(AbsTol))
    AbsTol=0;
end

if(nargin<9||isempty(RelTol))
    RelTol=0;
end

if(nargin<8||isempty(theta))
    theta=1;
end

if(nargin<7||isempty(maxIter))
    maxIter=2;
end

if(nargin<6||isempty(useNewton))
    useNewton=false;
end

if(nargin<5||isempty(givenVals))
    if(useNewton)
        [fCur,df]=f(x,curT);
    else
        fCur=f(x,curT);
    end
else
    fCur=[];
    df=[];
    if(isfield(givenVals,'fCur'))
        fCur=givenVals.fCur;
    end
    
    if(isfield(givenVals,'df'))
        df=givenVals.df;
    end
    
    if(isempty(fCur)||useNewton&&isempty(df))
        warning('givenVals provided but did not contain sufficient parameters. Calling f for the values.')
        if(useNewton)
            [fCur,df]=f(x,curT);
        else
            fCur=f(x,curT);
        end
    end
end

%Forward Euler step to get an initial estimate.
xCur=x+deltaT*fCur;
tNext=curT+deltaT;
if(maxIter==0&&nargout>1)
    didConverge=[];
    return
end

f0=fCur;
if(useNewton==false)
    %Use fixed point iteration.
    if(RelTol==0&&AbsTol==0)
        %If it should just iterate for the maximum time.
        for curIter=1:maxIter
            fCur=f(xCur,tNext);
            xCur=x+deltaT*(theta*fCur+(1-theta)*f0);
        end
        didConverge=[];
    else
        xOld=xCur;
        didConverge=false;
        for curIter=1:maxIter
            fCur=f(xCur,tNext);
            xCur=x+deltaT*(theta*fCur+(1-theta)*f0);
            
            diff=abs(xOld-xCur);
            if(all(diff<=AbsTol)||all(diff<=RelTol*abs(xCur)))
                didConverge=true;
                break; 
            end
            
            xOld=xCur;
        end
    end
else
    %Use a Newton iteration.
    xDim=size(x,1);
    I=eye(xDim,xDim);
    if(RelTol==0&&AbsTol==0)
        %If it should just iterate for the maximum time.
        %We already have the first fCur and df, so the first iteration is
        %done without calling f.
        F=xCur-x-deltaT*(theta*fCur+(1-theta)*f0);
        dF=I-theta*deltaT*df;
        xCur=xCur-dF\F;
        
        for curIter=2:maxIter
            [fCur,df]=f(xCur,tNext);
            F=xCur-x-deltaT*(theta*fCur+(1-theta)*f0);
            dF=I-theta*deltaT*df;
            xCur=xCur-dF\F;
        end
        didConverge=[];
    else
        xOld=xCur;
        didConverge=false;
        %We already have the first fCur and df, so the first iteration is
        %done without calling f.
        F=xCur-x-deltaT*(theta*fCur+(1-theta)*f0);
        dF=I-theta*deltaT*df;
        xCur=xCur-dF\F;
        diff=abs(xOld-xCur);
        if(all(diff<=AbsTol)||all(diff<=RelTol*abs(xCur)))
            didConverge=true;
            return; 
        end

        xOld=xCur;
        for curIter=2:maxIter
            [fCur,df]=f(xCur,tNext);
            F=xCur-x-deltaT*(theta*fCur+(1-theta)*f0);
            dF=I-theta*deltaT*df;
            xCur=xCur-dF\F;
            
            diff=abs(xOld-xCur);
            if(all(diff<=AbsTol)||all(diff<=RelTol*abs(xCur)))
                didConverge=true;
                return; 
            end
            
            xOld=xCur;
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
