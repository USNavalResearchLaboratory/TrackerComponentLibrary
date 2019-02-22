function Q=QTaylor(deltaT,xCur,curT,D,dadx,method)
%%QTAYLOR Get the process noise covariance under a nonlinear continuous-
%         time random process specified by the Langevin equation forward in
%         time by a step-size of deltaT using a Taylor scheme with additive
%         noise.
%
%INPUTS: deltaT The size of the single step over which to generate the
%               process noise covariance matrix.
%          xCur The initial target state at time curT.
%          curT The time of the initial state xCur.
%             D The diffusion function in the continuous-time stochastic
%               dynamic model. It takes the state and a time variable as
%               its arguments, D(x,t). Since the process noise must be
%               additive, D(x,t) should not depend on x.
%          dadx A xDimXxDim matrix of the derivative of the drift function
%               with respect to the state at state xCur and time curT.
%               This can also be a function that takes the state and a time
%               variable as its arguments, dadx(x,t).
%        method Set to 0 for the shorter Euler-Maruyama expansion, 1 for
%               the order 1.5 strong Taylor (default), and 2 for the order
%               2.0 weak Taylor.
%
%OUTPUTS: Q The process noise covariance matrix under a nonlinear
%           continuous-time random process specified by the Langevin
%           equation forward in time by a step-size of deltaT using a
%           strong Taylor scheme with additive noise.
%
%The process noise covariance matrix is derived from the stochastic order
%1.5 Taylor scheme with additive noise described in 10.4 of [1].
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isa(dadx,'function_handle'))
    dadx=dadx(xCur,curT);
end
if(nargin<6||isempty(method))
    method=1;
end

DCur=D(xCur,curT);

if(method==0)
    Q=deltaT*(DCur*DCur.');
else
    term1=deltaT*(DCur*DCur.');
    term2=(deltaT^2/2)*(DCur*(dadx*DCur).'+(DCur*(dadx*DCur).').');
    if(method==1)
        term3=(deltaT^3/3)*(dadx*DCur)*((dadx*DCur).');
    elseif(method==2)
        term3=(deltaT^3/4)*(dadx*DCur)*((dadx*DCur).');
    end
    
    Q=term1+term2+term3;
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
