function [c,E,exitVal]=minimaxPolyFit(f,xSpan,n,RelTol,maxIter)
%MINIMAXPOLYFIT Find an interpolating polynomial of order n to minimize
%               the maximum error approximating a scalar function f in the
%               range [xSpan(1);xSpan(2)]. The Remez algorithm is used.
%                
%INPUTS: f A handle to the function that is to be approximated. The use
%          y=f(x) take a scalar x and returns a value y.
%    xSpan The range of inputs of x over which a minimax interpolating
%          polynomial is desired. The format is [xMin;xMax]. If this is
%          omitted or an empty matrix is passed, the values [-1;1] are
%          used.
%        n The degree of the polynomial approximation desired. If this is
%          omitted or an empty matrix is passed, the default value of n=5
%          is used. If too small a value is used and the function
%          oscillates a lot, the output polynomial will not be a true
%          minimax polynomial, because the algorithm assumes that there is
%          one maxima between regions chosen for selecting points. Also,
%          if f is an mth order polynomial, but is symmetric, then n
%          should be m+1 due to assumptions in the algorithm (the
%          resulting output will have the extra coefficient zero).
%   RelTol The relative tolerance of the maximum error used for
%          determining convergence. Iterations adjust the maximum error
%          observed with respect to points at which the function is
%          evaluated. When the maximum error stops changing across
%          iterations, convergence is assumed. The relative maximum error
%          is abs(E-EPrev)/abs(E); where E is the error and EPRev is the
%          error from the previous iteration. If this parameter is omitted
%          or an empty matrix is passed, a value of 1e-12 is used.
%  maxIter The maximum number of iterations for convergence. If omitted, a
%          value of 50 is used.
%
%OUTPUTS: c An interpolating polynomial such that polyval(c,xRel) gives the
%           interpolated value of f where xRel is the fractional distance
%           between xSpan(1) and xSpan(2). The value xRel ranges from 0 to
%           1. The development of the interpolation across 0 and 1 is often
%           numerically better than explicitly making the polynomial for
%           values from xSpan(1) to xSpan(2).
%         E The maximum error of the interpolating polynomial.
%   exitVal A value indicating how the algorithm terminated. Possible
%           values are:
%           0: The tolerance goal was met.
%           1: the maximum number of iterations was met.
%
%The algorithm used to fit a polynomial to an arbitrary nonlinear function
%is Remez algorithm. This is described ine Remez algorithm is well known
%and is described in Chapter 6 of [1].
%
%To initalize the Remez algorithm, one must choose a set of points in the
%region to be interpolated. One might be tempted to unifromely space the
%points, but as noted in [2] and Chaper 8.3 of [3], that can be a very bad
%idea. Rather, the use of Chebychev nodes is better and has much better
%asymptotic properties as the number of nodes increases. 
%
%Given a set of points in xSpan, the algorithm find an interpolating
%polynomial that mades the magnitude of the error at all of the points the
%same and the sign of the error alternate. This is because there is a
%theorem that such an alternating equierror solution is a minimax solution.
%However, only a minimax solution for the initial set of points is given.
%Thus, the algorithm tries to move the points to maximize the error.
%
%Since the sign of the error between any two points flips, there is at
%least one zero between the points. This algorithm assumes that there is
%exactly one zero between the points, which should be reasonable if n is
%chosen to be large enough. In such an instance, the true error maxima of
%the function are between the zeros (and the outer bounds) of the function.
%Thus, the points can be moved there. Basically by maximizing the error at
%the chosen points, and getting an interpolating polynomial with equal
%error at the points, one finds the minimax approximation.
%
%As an example, consider interpolating the exponential function between 0
%and 2 with a sixth-order polynomial:
% f=@(x)exp(x);
% xSpan=[0;2];
% n=6;
% [c,E]=minimaxPolyFit(f,xSpan,n)
% %The maximum error E should be about 8.7e-6.
% %The resulting interpolating polynomial can be plotted as a fraction of
% %the distance between the points in xSpan. We will plot the error.
% numPoints=500;
% x=linspace(xSpan(1),xSpan(2),numPoints);
% xFrac=(x-xSpan(1))/(xSpan(2)-xSpan(1));
% error=polyval(c,xFrac)-exp(x);
% plot(x,error)
%
%REFERENCES:
%[1] V. K. Dzyadyk and I. A. Shevchuk, Theory of Uniform Approximation
%    of Functions by Polynomials, 1st ed. Berlin: Walter de Gruyter, 2008.
%[2] R. B. Platte, L. N. Trefethen, and A. B. Kuijlaars, "Impossibility of
%    fast stable approximation of analytic functions from equispaced
%    samples," SIAM Review, vol. 53, no. 2, pp. 308-318, 2011.
%[3] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Australia:
%    Brooks/ Cole, 2011.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(maxIter))
   maxIter=50; 
end

if(nargin<4||isempty(RelTol))
   RelTol=1e-12; 
end

if(nargin<3||isempty(n))
    n=5;
end

if(nargin<2||isempty(xSpan))
    xSpan=[-1;1];
end

%Step 1: Choose an initial set of nodes for interpolation. 
%n+2 nodes are needed for an order-n polynomial approximation. The initial
%guess is set to the roots of the Chebyshev polynomial of the first kind.
%Such a choice is standard in most implementations of the Remez algorithm,
%because the Chebyshev polynomials distribute more nodes to the outside of
%the interval then to the center, reducing Runge's phenomenon. 
k=(n+2):-1:1;
%The Chebyshev nodes for the (-1,1) span.
xNodes=cos((2*k(:)-1)/(2*(n+2))*pi);
%We are reducing the span of the nodes from -1/2 to 1/2 so that the
%interpolating polynomial on the output of the function takes an input from
%0 to 1.
xNodes=xNodes/2;

curIter=1;
EPrev=Inf;%This is used to determine convergence of the error estimate E.
while(1)
    %Step 2: Solve for the maximum error E and the coefficients of the
    %        interpolating polynomial P given the current set of nodes.
    %        This involves solving a linear system of equations.

    %Allocate space for the matrix.
    hMat=zeros(n+2,n+2);

    xDelta=(xNodes+1/2);
    hMat(:,1)=1;%The zeroth-order terms.
    for curCol=2:(n+1)
        hMat(:,curCol)=hMat(:,curCol-1).*xDelta;%This is xDelta.^(curCol-1)
    end

    %The last column corresponds to the E error term. This last column is
    %why the matrix is not a Vandermonde matrix.
    hMat(:,n+2)=(-1).^(0:(n+1));

    %Map the (-1,1) Chebyshev nodes into the range given by xSpan.
    xFNodes=(xSpan(2)-xSpan(1))*(xNodes+1/2)+xSpan(1);
    
    %Solve the system for the polynomial coefficients c and the error E
    coeffs=hMat\f(xFNodes);

    %The flipud is to make the ordering the same as that used by the
    %polyval function.
    c=flipud(coeffs(1:(n+1)));
    E=coeffs(n+2);
    
    relTolVal=abs(E-EPrev)/abs(E);
    
    %Check for convergence. Convergence is determined is the E value stops
    %changing significantly.
    %The check for the relative tolerance value being a nan deals with the
    %case where E=0 and E-EPrev=0.
    if(isnan(relTolVal)||relTolVal<RelTol)
        E=abs(E);
        exitVal=0;
        break; 
    end
    
    %Return if the maximum number of iterations has passed.
    if(curIter==maxIter)
        E=abs(E);
        exitVal=1;
        break;
    end
    
    EPrev=E;

    %The coefficients c form the order-n minimax approximation of the 
    %function at the n nodes. The error at those nodes has magnitude E.
    %However, these are most likely not the maximum error points for the
    %full function f. We want to move the points xNodes so that we are
    %minimizing the maximum error with respect to f in the interval.

    %Step 3: The Exchange Step. Given the current polynomial interpolation,
    %move the points to where it maximizes the magnitude of the error
    %function.

    %Step 3.1: Find the zeros of the error function between the control
    %          points.
    costFunZero=@(x)errorFun(x,xSpan,c,f);
    %There are n+1 points; the extra 2 points allocated are to make finding
    %the extrema in the next part of the step simpler. 
    xIntervals=zeros(n+3,1);
    for curPoint=1:(n+1)
       %Find the zero of the error function between points curPoint and
       %curPoint+1.
       %Precision is arbitrary, but since X is from -1/2 to 1/2, it makes
       %sense to hard-code in the tolerance on X (and use a higher
       %tolernace than is the default)
       options=optimset('TolX',1e-14);
       xIntervals(curPoint+1)=fzero(costFunZero,[xNodes(curPoint),xNodes(curPoint+1)],options);
    end

    %Step 3.2: Find the extreme points and make them the new nodes.
    %The extrema are bracketed between each pair of roots. There are also
    %two more extrema located between the start of the interval xSpan(1)
    %and the first root and between the last root and the end of the
    %interval xSpan(2). We want to find all of the extrema. We will do line
    %searches to find those points. Those points are placed into xNodes as
    %they will form the new nodes for the next iteration.
    xIntervals(1)=-1/2;
    xIntervals(end)=1/2;

    costFunMax=@(x)negAbsErrorFun(x,xSpan,c,f);
    for curPoint=1:(n+2)
        %Find the maximum of the error function in the bounded interval 
        %using a line search algorithm.
        options=optimset('TolX',1e-14);
        xNodes(curPoint)=fminbnd(costFunMax,xIntervals(curPoint),xIntervals(curPoint+1),options);
    end
    
    %Step 4: Iterate.
    curIter=curIter+1;
end
end

function val=errorFun(x,xSpan,c,f)
   %Map the (-1/2,1/2) Chebyshev node into the range given by xSpan.
   xf=(xSpan(2)-xSpan(1))*(x+1/2)+xSpan(1);

   val=polyval(c,x+1/2)-f(xf);
end

function val=negAbsErrorFun(x,xSpan,c,f)
   %Map the (-1/2,1/2) Chebyshev node into the range given by xSpan.
   xf=(xSpan(2)-xSpan(1))*(x+1/2)+xSpan(1);

   val=-abs(polyval(c,x+1/2)-f(xf));
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
