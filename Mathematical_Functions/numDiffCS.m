function J=numDiffCS(x,f,fDim,deltaX)
%%NUMDIFFCS Perform numeric differentiation of a real function of a real
%           variable when the function happens to be complex analytic by
%           using the complex-step approximation. This can be dramatically
%           more accurate than the finite-difference approximations used in
%           the numDiff function. The technique is derived in [1].
%
%INPUTS: x The xDimX1 real vector or scalar point at which the derivative
%          of the (possibly vector) function is desired.
%        f The scalar or vector function that is to be differentiated. The
%          function f must take x as its parameter and its output should be
%          a scalar or a column vector. The function should be complex
%          analytic but return a real value when given a real input.
%     fDim The dimensionality of the output of f. If this is omitted, or
%          an empty matrix is passed, it is assume dthat fDim=1.
%   deltaX Optionally, the step size can be specified. This is either a
%          scalar or an xDimX1 vector. The default if omitted or an empty
%          matrix is passed is 1e-50.
%
%OUTPUTS: J An fDimXxDim Jacobian matrix. Each column is the derivative
%           vector of f with respect to the corresponding element of x. If
%           at any point the function f returned an empty matrix, then J
%           will be an empty matrix.
%
%EXAMPLE:
%This example makes f a vector function using the two example functions
%mentioned in [1]. Since we can solve for the derivative analytically, we
%compute it directly from the solution using double precision (the Matlab
%default). We then compare that to a high precision solution that was
%obtained using extended precision arithmetic. We also find the error in
%using numDiff compared to numDiffCSErr. The result from numDiffCSerr is
%about as accurate as direct computation, and is multiple orders of
%magnitude more accurate than using numDiff.
% f=@(x)[x^(9/2);
%        exp(x)/(sin(x)^3+cos(x)^3)];
% x0=1.5;
% df=@(x)[(9/2)*x^(7/2);
%          exp(x)*(2*cos(3*x)+3*sin(x)+sin(3*x))/(2*(cos(x)^3+sin(x)^3)^2)];
% derivValNumeric=df(x0);
% derivValHighPrec=[18.6008127342597586831856259422978627577414136812365887877858;
%                   3.62203370071632604259777098114019696227123307201775463275874];
% numDiffVal=numDiff(x0,f,2);
% complexDiffVal=numDiffCS(x0,f,2);
% directError=norm(derivValNumeric-derivValHighPrec)
% numDiffError=norm(numDiffVal-derivValHighPrec)
% numDiffCSError=norm(complexDiffVal-derivValHighPrec)
%
%REFERENCES:
%[1] W. Squire and G. Trapp, "Using complex variables to estimate
%    derivatives of real functions," SIAM Review, vol. 40, no. 1,
%    pp. 110-112, Mar. 1998.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(fDim))
    fDim=1;
end

xDim=length(x);

if(nargin<4||isempty(deltaX))
    deltaX=1e-50*ones(xDim,1);
end

J=zeros(fDim,xDim);
for curDim=1:xDim
    x(curDim)=x(curDim)+1i*deltaX(curDim);
    val=imag(f(x)./deltaX(curDim));
    if(isempty(val))
       J=[];
       return;
    end
    J(:,curDim)=val;
    x(curDim)=real(x(curDim));
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
