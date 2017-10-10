function [a,t]=BSplinePolyFit(tau,y,k,numCentralKnots)
%%BSPLINEPOLYFIT Given a set of real scalar points and real or complex
%        scalar function values at those points, this function returns the
%        knots and coefficients for a set of b-splines that can be used to
%        interpolate to a desired polynomial degree within the sampled
%        region of the function. A least squared interpolant fit is also
%        possible. B-splines are useful, because they can have continuous
%        derivatives across the entire sampled region and they can be
%        easily differentiated. This function uses the heuristic
%        approximation for the approximate knots from Equation 10 of
%        Chapter 13 of [1] rather than finding the optimal knots using
%        Newton's method. The first and last knots are repeated and
%        coincide with the tau(1) and tau(end). This function can provide
%        interpolation coefficients for multiple sets of values at once. 
%
%INPUTS: tau A numPointsX1 or 1 1XnumPoints vector of points where
%          the function is sampled. These must be provided sorted in
%          increasing order and must be real.
%        y The numPointsXnumSets sets of function values at the points in
%          tau that are to be matched. Each set produces a corresponding
%          set of interpolation coefficients. The sets are not related. the
%          values in y can be complex.
%        k The order of the b-spline. The value k-1 is the polynomial order
%          of the approximation being performed. If this parameter is
%          omitted or an empty matrix is passed, the default of k=4 is
%          used. This value must be real.
% numCentralKnots This optional parameter is used if one passes a higher
%          density of points in tau than is needed, because one wants to
%          perform least squares interpolation. This is the number of knots
%          in the central region of knots. It varies from 0 to n-k. The
%          higher the number, the more knots are used total. A total of
%          N=numCentralKnots+2*k will be returned by this function, though,
%          as explained below, the first and last knots are not used.
%
%OUTPUTS: a A numAPointsXnumSets collections of coefficients for the
%           b-splines covering all of the interpolation intervals for each
%           set of values to be matched. This and the following output can
%           be passed to BSplineInterpVal to perform interpolation over the
%           region tau(1) to tau(end). if y is complex, then a is complex.
%         t The NX1 collection of knots for the interpolation function. The
%           first and last k-1 knots are outside of the ends of the region
%           with data or mark the ends of the region with data. These
%           values are the same for all sets of data that are interpolated
%           in a. These values are real.
%
%Equation 7 in Chapter IX of [1] expresses an interpolated function in
%terms of a weigted combination of b-spline polynomials. Solving for the
%weights is simply a matter of solving a set of linear equations. However,
%constructing the set of linear equations requires choosing appropriate
%knots. Here, we use the approximate knots from Equation 10 of Chapter 13
%of [1] with the first k knots set to tau(1) and the last k knots set to
%tau(numPoints).
%
%The ability to do a least-squares fit involves making fewer inner knots
%than suggested in the text. Rather than making n-k knots, we make
%numCentralKnots. However, this invalidates the approximation in Equation
%10 of Chapter 13 of [1]. The correction applied is to make the average
%taken in the equation to involve more points.
%
%EXAMPLE 1:
%Here, we fit a curve to a fourth-order function using a series of third
%-order functions. The fit is mediocre. The irregular sample interval
%contributes to the error.
% f=@(t)(t.^4-2*t.^2+t);
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=f(tau);
% k=4;%k-1=3 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500);
% yp=BSplineInterpVal(x,t,a,k);
% y=f(x);
% figure(1)
% clf
% hold on
% plot(x,y,'-k','linewidth',4)
% plot(x,yp,'--r','linewidth',2)
%
%EXAMPLE 2:
%This is the same as example 1, except we use fourth-order polynomial
%segments to interpolate the fourth-order function. Now, the fit is
%perfect.
% f=@(t)(t.^4-2*t.^2+t);
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=f(tau);
% k=5;%k-1=4 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500);
% yp=BSplineInterpVal(x,t,a,k);
% y=f(x);
% figure(2)
% clf
% hold on
% plot(x,y,'-k','linewidth',4)
% plot(x,yp,'--r','linewidth',2)
%
%EXAMPLE 3:
%In this example, a least-squares fit is used for the interpolation. One
%will see that the least squares fitted curve (red) looks much better than
%the exactly fitted curve (blue) when compared to the truth (black).
% f=@(t)(t.^4-2*t.^2+t);
% tau=-2:0.25:2;
% y=f(tau)+0.5*randn(size(tau));
% y=y(:);
% k=4;%k-1=3 is the order....
% 
% %Least-squares fit on the noisy data.
% numCentralKnots=2;
% [aLS,tLS]=BSplinePolyFit(tau,y,k,numCentralKnots);
% 
% %Standard fit on the noisy data.
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500);
% yp=BSplineInterpVal(x,t,a,k);
% yLS=BSplineInterpVal(x,tLS,aLS,k);
% y=f(x);
% figure(2)
% clf
% hold on
% plot(x,y,'-k','linewidth',4)
% plot(x,yLS,'--r','linewidth',2)
% plot(x,yp,'--c','linewidth',2)
%
%EXAMPLE 4:
%In this example two sets of values for interpolating different functions
%are given at once and then the fits to the two functions are both plotted.
%Interpolating sets of functions is useful when one wishes to interpolate a
%function that has a multidimensional output.
% f1=@(t)(t.^4-2*t.^2+t);
% f2=@(t)(t.^3-2*t.^2+t);
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=[f1(tau),f2(tau)];
% k=4;%k-1=3 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500).';
% yp=BSplineInterpVal(x,t,a,k);
% y=[f1(x),f2(x)];
% %Plot the first interpolated function
% figure(1)
% clf
% hold on
% plot(x,y(:,1),'-k','linewidth',4)
% plot(x,yp(:,1),'--r','linewidth',2)
% %Plot the second interpolated function
% figure(2)
% clf
% hold on
% plot(x,y(:,2),'-k','linewidth',4)
% plot(x,yp(:,2),'--r','linewidth',2)
%
%EXAMPLE 5:
%In this example, we fit a complex function. We then look at how well each
%component fits.
% numPoints=500;
% xMin=0;
% xMax=10;
% x=linspace(xMin,xMax,numPoints);
% f=@(x)exp(-1j*2*x).*(4*x.^3+x.*exp(1j*2*x.^2)+12*x+x.^3.*exp(-1j*2*x.^2));
% 
% fx=f(x);
% k=5;
% [a,t]=BSplinePolyFit(x,fx(:),k);
% 
% %Interpolate on a finer grid.
% numPoints=2000;
% x=linspace(xMin,xMax,numPoints);
% fx=f(x);
% fxInterp=BSplineInterpVal(x,t,a,k);
% 
% %Plot the function values and the interpolated values for the real and
% %complex parts.
% figure(1)
% clf
% hold on
% plot(x,real(fx),'-k','linewidth',2)
% plot(x,imag(fx),'--k','linewidth',2)
% plot(x,real(fxInterp),'-r')
% plot(x,imag(fxInterp),'--r')
% 
% %Plot the interpolation error.
% figure(2)
% clf
% hold on 
% plot(x,real(fxInterp(:))-real(fx(:)))
% plot(x,imag(fxInterp(:))-imag(fx(:)))
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(tau);

if(nargin<3||isempty(k))
   k=4; 
end

if(nargin<4||isempty(numCentralKnots))
    numCentralKnots=n-k;
elseif(numCentralKnots<0||numCentralKnots>n-k)
    error('The value of numCentralKnots is invalid.')
end

deltaVal=n-k-numCentralKnots;

t=zeros(numCentralKnots+2*k,1);
t(1:(k))=tau(1);
t((numCentralKnots+k+1):(numCentralKnots+2*k))=tau(n);

%The ad-hoc point choice from 10 of Chapter 13, modified to take more
%points if a least-squares solution is being performed.
for i=1:numCentralKnots
    t(k+i)=sum(tau((i+1):(i+k-1+deltaVal)))/(k-1+deltaVal);
end

%Given the locations of the knots, we build the equations 
B=zeros(n,numCentralKnots+k);

span=1:k;
tIdx=k;
for m=1:n
    if(tau(m)>t(tIdx+1))
        tIdx=tIdx+1;
        span=span+1;
    end

    B(m,span)=evalBSplinePolys(t,tau(m),k,tIdx);
end

a=B\y;

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
