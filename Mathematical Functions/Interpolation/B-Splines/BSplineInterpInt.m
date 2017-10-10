function [t,a,k]=BSplineInterpInt(t,a,k)
%%BSPLINEINTERPINT Given a set of b-spline interpolation weights a and
%            knots t for scalar interpolation starting at the point x0 assu
%            assuming that t(1) is repeated k times,(e.g. the smallest tau
%            value used in the BSplinePolyFit function was x0), compute the
%            weights and knots to interpolate the definite integral of the
%            curve from x0 to a given point using b-spline interpolation.
%            If t(1) is not repeated, then this can still be used, but one
%            should consider it an indefinite integral with a particular
%            additive constant and should difference two of these to
%            perform a definite integral.
%
%INPUTS: t The NX1 or 1XN set of knots for the interpolation function. The
%          first and last k-1 knots are outside of the ends of the region
%          with data or mark the ends of the region with data. This must be
%          real.
%        a The npXnumSets collection of b-spline coefficients for
%          interpolating over a certain region for numSets sets of
%          interpolation problems. Such coefficients could be obtained from
%          the function BSplinePolyCoeffs along with t if one wishes to
%          interpolate over gridded data. This can be real or complex.
%        k The order of the B-spline polynomials. The order of the
%          interpolation accuracy is k-1. 
%
%OUTPUTS: a The (np+1)XnumSets modified coefficients.
%         t The (N+2)X1 modified knots.
%         k The modified B-spline order.
%
%This function implements Equation 22 in Chapter X of [1]. The modification
%to t simply reflects the modification to the overall length of a.
%
%EXAMPLE 1:
%Here, we plot the definite integral of the derivative a function versus
%the b--interpolated value. In this example, we chose the order to be high
%enough for the fit to be exact.
% f=@(t)(t.^4-2*t.^2+t);
% %The definite integral will be from -2
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=f(tau);
% k=5;%k-1=4 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% [t,a,k]=BSplineInterpInt(t,a,k);
% numPoints=100;
% x=linspace(-2,2,numPoints);
% y=-(14/15)+x.^2/2-(2*x.^3)/3+x.^5/5;
% 
% yInterp=BSplineInterpVal(x,t,a,k);
% figure(1)
% clf
% hold on
% plot(x,y,'-k','linewidth',4);
% plot(x,yInterp,'--r','linewidth',2);
%
%EXAMPLE 2:
%This is similar to the previous example, except multiple sets of
%coefficients are computed and evaluated at once.
% f1=@(t)(t.^4-2*t.^2+t);
% f2=@(t)(t.^3-2*t.^2+t);
% %The definite integral will be from -2
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=[f1(tau),f2(tau)];
% k=5;%k-1=4 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% [t,a,k]=BSplineInterpInt(t,a,k);
% numPoints=100;
% x=linspace(-2,2,numPoints).';
% y=[-(14/15)+x.^2/2-(2*x.^3)/3+x.^5/5,-(34/3)+x.^2/2-(2*x.^3)/3+x.^4/4];
% 
% yInterp=BSplineInterpVal(x,t,a,k);
% figure(1)
% clf
% hold on
% plot(x,y(:,1),'-k','linewidth',4);
% plot(x,yInterp(:,1),'--r','linewidth',2);
% 
% figure(2)
% clf
% hold on
% plot(x,y(:,2),'-k','linewidth',4);
% plot(x,yInterp(:,2),'--r','linewidth',2);
%
%EXAMPLE 3:
%This is an example of where a complex function parameterized by a real
%value is interpolated.
% numPoints=500;
% xMin=0;
% xMax=3;
% x=linspace(xMin,xMax,numPoints);
% f=@(x)(4*x.^3+x.*exp(1j*2*x.^2)+12*x+x.^3.*exp(-1j*2*x.^2)+10*1i*cos(x));
% 
% %The integral is from 0 to x.
% fIntx=@(x)(-(1/4)*1j*(-1+exp(2*1j*x.^2))+6*x.^2+x.^4+(1/8)*(-1+exp(-2*1j*x.^2).*(1+2*1j*x.^2))+10*1i.*sin(x));
% 
% fx=f(x);
% k=5;
% [a,t]=BSplinePolyFit(x,fx(:),k);
% [t,a,k]=BSplineInterpInt(t,a,k);
% 
% %Interpolate on a finer grid.
% numPoints=2000;
% x=linspace(xMin,xMax,numPoints);
% fx=fIntx(x);
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

N=length(t);
numA=size(a,1);
numSets=size(a,2);

betaVal=zeros(numA+1,numSets);
for i=1:numA
    betaVal(i+1,:)=sum(bsxfun(@times,a(1:i,:),(t((1:i)+k)-t(1:i))),1)/k;
end

a=betaVal;
k=k+1;
t=[t(1);t(:);t(N)];

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
