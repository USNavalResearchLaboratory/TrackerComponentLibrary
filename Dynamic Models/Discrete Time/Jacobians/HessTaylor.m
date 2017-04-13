function F=HessTaylor(deltaT,dadx,d2adx,method)
%%HESSTAYLOR Simulate a nonlinear continuous-time random process specified
%            by the Langevin equation forward in time by a step-size of
%            deltaT using a Taylor scheme.
%
%INPUTS: deltaT The size of the single step over which to generate the
%               state transition matrix.
%          dadx A xDimXxDim matrix of the derivative of the drift function
%               with respect to the state at state xCur and time curT.
%         d2adx A xDimXxDimXxDim matrix of the second derivative of the
%               drift function with respect to the state at state xCur and
%               time curT. The value at point (m,k,l) represents
%               d2a(m)/dx(k)dx(l). If not provided, this is assumed to be
%               zero.
%        method Set to 0 for the shorter Euler-Maruyama expansion,
%               otherwise will the Taylor expansion will be used (default).
%
%OUTPUTS: H The Hessian matrix of a nonlinear continuous-time random
%           process specified by the Langevin equation forward in time by a
%           step-size of deltaT using a strong Taylor scheme.
%
%The second derivative state prediction matrix is derived from the
%stochastic order 1.5 Taylor scheme described in 10.4 of [1].
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(dadx,1);
if(nargin<4||isempty(method))
    method=1;
end

if(method==0)
    F=deltaT*d2adx;
else
    %Assumes 3rd derivative is zero
    term1=zeros(xDim,xDim,xDim);
    term2=zeros(xDim,xDim,xDim);
    term3=zeros(xDim,xDim,xDim);
    for m=1:xDim
        term1(:,:,m)=d2adx(:,:,m)*dadx;
        term2(:,m,:)=d2adx(:,:,m)*dadx;
        term3(:,:,m)=dadx*d2adx(:,:,m);
    end
    
    F=deltaT*d2adx+(deltaT^2/2)*(term1+term2+term3);
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
