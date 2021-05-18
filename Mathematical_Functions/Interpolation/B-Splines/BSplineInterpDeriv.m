function [t,a,k]=BSplineInterpDeriv(t,a,k,numDerivs)
%%BSPLINEINTERPDERIV Given a set of b-spline interpolation weights a and
%            knots t for scalar interpolation, compute the weights and
%            knots to interpolate the numDerivs derivative value using
%            b-splines interpolation.
%
%INPUTS: t The NX1 or 1XN set of knots for the interpolation function. The
%          first and last k-1 knots are outside of the ends of the region
%          with data or mark the ends of the region with data. These vslues
%          must be real.
%        a The npXnumSets collection of b-spline coefficients for
%          interpolating over a certain region for numSets sets of
%          interpolation problems. Such coefficients could be obtained from
%          the function BSplinePolyCoeffs along with t if one wishes to
%          interpolate over gridded data. These values can be complex.
%        k The order of the B-spline polynomials. The order of the
%          interpolation accuracy is k-numDerivs-1.
% numDerivs The number of derivatives to take. This must be less than k-1.
%
%OUTPUTS: a The (np-numDerivs)XnumSets modified coefficients.
%         t The (N-2*numDerivs)X1 modified knots.
%         k The modified B-spline order.
%
%This function implements Equation 12b in Chapter X of [1]. The
%modifications to t simply reflect the modification to the overall length
%of a.
%
%The first and last elements in t are not used in this function. However,
%the equations in [1] assume that those elements are present as they end up
%being used when computing integrals.
%
%EXAMPLE 1:
%Here, we approximate the derivative of a fourth-order function using
%piecewise third-order b-splines and plot the results.
% f=@(t)(t.^4-2*t.^2+t);
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=f(tau);
% k=4;%k-1=3 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500);
% numDerivs=1;
% [t,a,k]=BSplineInterpDeriv(t,a,k,numDerivs);
% ypDeriv=BSplineInterpVal(x,t,a,k);
% yDeriv=1-4*x+4*x.^3;%First derivative
% 
% figure(1)
% clf
% hold on
% plot(x,yDeriv,'-k','linewidth',4)
% plot(x,ypDeriv,'--r','linewidth',2)
%
%EXAMPLE 2:
%This is similar to the previous example, except multiple sets of
%coefficients are computed and evaluated at once.
% f1=@(t)(t.^4-2*t.^2+t);
% f2=@(t)(t.^3-2*t.^2+t);
% tau=[-2;-1;0;0.5;1;1.5;2];
% y=[f1(tau),f2(tau)];
% k=4;%k-1=3 is the order....
% [a,t]=BSplinePolyFit(tau,y,k);
% 
% x=linspace(-2,2,500).';
% numDerivs=1;
% [t,a,k]=BSplineInterpDeriv(t,a,k,numDerivs);
% ypDeriv=BSplineInterpVal(x,t,a,k);
% yDeriv=[1-4*x+4*x.^3,1-4*x+3*x.^2];%First derivative
% 
% figure(1)
% clf
% hold on
% plot(x,yDeriv(:,1),'-k','linewidth',4)
% plot(x,ypDeriv(:,1),'--r','linewidth',2)
% 
% figure(2)
% clf
% hold on
% plot(x,yDeriv(:,2),'-k','linewidth',4)
% plot(x,ypDeriv(:,2),'--r','linewidth',2)
%
%EXAMPLE 3:
%This is an example where a complex function parameterized by a single real
%value is differentiated and interpolated.
% numPoints=300;
% xMin=0;
% xMax=5;
% x=linspace(xMin,xMax,numPoints);
% f=@(x)exp(-1j*2*x).*(4*x.^3+x.*exp(1j*2*x.^2)+12*x+x.^3.*exp(-1j*2*x.^2));
% 
% fDx=@(x)exp(-2*1j*x.*(1+x)).*(exp(4*1j*x.^2).*(1+2*1j*x.*(-1+2*x))+x.^2.*(3-2*1j*x.*(1+2*x))+4*exp(2*1j*x.^2).*(3+x.*(-6*1j+(3-2*1j*x).*x)));
% 
% fx=f(x);
% k=5;
% [a,t]=BSplinePolyFit(x,fx(:),k);
% numDerivs=1;
% [t,a,k]=BSplineInterpDeriv(t,a,k,numDerivs);
% 
% %Interpolate on a finer grid.
% numPoints=2000;
% x=linspace(xMin,xMax,numPoints);
% fx=fDx(x);
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

numA=size(a,1);

for curDeriv=1:numDerivs
    aNext=a;
    for idxA=(curDeriv+1):numA
        denom=(t(idxA+k-curDeriv)-t(idxA))/(k-curDeriv);

        if(denom==0)
            aNext(idxA,:)=0;
        else
            aNext(idxA,:)=(a(idxA,:)-a(idxA-1,:))/denom;
        end
    end

    a=aNext;
end
a=a((1+numDerivs):end,:);
t=t((1+numDerivs):(end-numDerivs));
k=k-numDerivs;
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
