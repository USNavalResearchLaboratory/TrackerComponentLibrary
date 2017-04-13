function F=FTaylor(deltaT,xCur,curT,a,dadx,d2adx,method)
%%FTAYLOR Simulate a nonlinear continuous-time random process specified by
%         the Langevin equation forward in time by a step-size of deltaT
%         using a Taylor scheme.
%
%INPUTS: deltaT The size of the single step over which to generate the
%               state transition matrix.
%          xCur The initial target state at time curT.
%          curT The time of the initial state xCur.
%             a The drift function in the continuous-time stochastic
%               dynamic model. It takes the state and a time variable as
%               its arguments, a(x,t).
%          dadx A xDimXxDim matrix of the derivative of the drift function
%               with respect to the state at state xCur and time curT.
%               This can also be a function that takes the state and a time
%               variable as its arguments, dadx(x,t).
%         d2adx A xDimXxDimXxDim matrix of the second derivative of the
%               drift function with respect to the state at state xCur and
%               time curT. The value at point (m,k,l) represents
%               d2a(m)/dx(k)dx(l). This can also be a function that takes
%               the state and a time variable as its arguments,
%               d2adx(x,t). If not provided, this is assumed to be zero.
%        method Set to 0 for the shorter Euler-Maruyama expansion, 
%               otherwise will use Taylor expansion (default).
%
%OUTPUTS: F The state transition matrix under a nonlinear continuous-time
%           random process specified by the Langevin equation forward in
%           time by a step-size of deltaT using a strong Taylor scheme.
%
%The state prediction matrix is derived from the stochastic order 1.5
%Taylor scheme described in 10.4 of [1].
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

aCur=a(xCur,curT);
xDim=size(aCur,1);

if(nargin<6 || isempty(d2adx))
    d2adx=zeros(xDim,xDim,xDim);
end
if(nargin<7||isempty(method))
    method=1;
end
if(isa(dadx,'function_handle'))
    dadx=dadx(xCur,curT);
end
if(isa(d2adx,'function_handle'))
    d2adx=d2adx(xCur,curT);
end

if(method==0)
    F=eye(xDim)+deltaT*dadx;
else
    %assumes 3rd derivative is zero and all cross-derivatives are zero
    term1=zeros(xDim);
    for n=1:xDim
        term1(:,n)=d2adx(:,:,n)*aCur;
    end
    F=eye(xDim)+deltaT*dadx+(deltaT^2/2)*(term1+dadx*dadx);
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
