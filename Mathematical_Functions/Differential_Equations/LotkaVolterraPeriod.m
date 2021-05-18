function [T,totalError,exitCode]=LotkaVolterraPeriod(x0,y0,a,b,c,d,intParams)
%%LOTKAVOLTERRAPERIOD Determine the period of a given set of Lotka-Volterra
%           equations. These are also known as predator-prey equations.
%           These are the following set of two coupled differential
%           equations:
%           dx=a*x-b*x*y
%           dy=c*x*v-d*y
%           with a,b,c,d>0 and x,y with initial conditions >0. Such
%           equations are useful for analyzing the performance of
%           techniques for solving initial values problems.
%
%INPUTS: x0, y0 The real initial conditions of the two coupled processes.
%               x0,y0 >0.
%       a,b,c,d The real coefficients of the model, all >0.
%     intParams Optionally, one can pass a structure affecting how the
%               numerical integration in this function is performed. This
%               can be a structure with fields RelTol, AbsTol, n, and
%               maxSearchReg having respective default values of 1e-8,
%               1e-11, 8, and 1000. These values correspond to the
%               same-named parameters of the integral1DAdaptive function.
%
%OUTPUTS: T The period of the Lotka-Volterra system.
% totalError An estimate of the total error of the solution. This is
%           derived from the same outpuf of the integral1DAdaptive
%           function.
%  exitCode The exit code returned by the integral1DAdaptive function. 0
%           means that the absolute ore relative error condition was
%           fulfilled and thus, the algorithm converged.
%
%This function implements the integral given by Theorem 6.2 of [1].
%
%EXAMPLE:
%Here, we determine the period of s system and then we perform numerical
%integration over that period to verify that the solution looks like a
%closed loop.
% a=2/3;
% b=4/3;
% c=1;
% d=1;
% x0=1/5;
% y0=1/5;
% [T,totalError,exitCode]=LotkaVolterraPeriod(x0,y0,a,b,c,d)
% %The Lotka-Volterra differential system.
% f=@(x,t)([a;-d].*x+[-b;c]*prod(x));
% numSteps=1e4;
% xStart=[x0;y0];
% theTimes=linspace(0,T,numSteps);
% xList=RungeKAtTimes(xStart,theTimes,f,Inf);
% figure(1)
% clf
% plot(xList(1,:),xList(2,:),'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('y');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the plot is a closedl oop and that reducing the
%integration time a little bit opens the loop, indicating that the solution
%was correct.
%
%REFERENCES:
%[1] S.-D. Shih, "The period of a Lotka-Volterra system," Taiwanese Journal
%    of Mathematics, vol. 1, no. 4, pp. 451-470, Dec. 1997.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

RelTol=1e-8;
AbsTol=1e-11;
n=8;
maxSearchReg=1000;

if(nargin>6&&~isempty(intParams))
    if(isfield(intParams,'RelTol'))
        RelTol=intParams.RelTol;
    end
    
    if(isfield(intParams,'AbsTol'))
        AbsTol=intParams.AbsTol;
    end
    
    if(isfield(intParams,'n'))
        n=intParams.n;
    end
    
    if(isfield(intParams,'maxSearchReg'))
        n=intParams.maxSearchReg;
    end
end

%The energy of the system, Equation 1.5. E>=0.
E=max(0,c*x0-d+b*y0-a-a*log((b/a)*y0)-d*log((c/d)*x0));

func=@(s)intArg(s,E,a,d);

[intVal,totalError,exitCode]=integral1DAdaptive(func,[0;E],n,2,[],RelTol,AbsTol,maxSearchReg);
T=(1/(a*d))*intVal;
totalError=(1/(a*d))*totalError;
end

function val=intArg(s,E,a,d)
%%INTARG The argument fo the integral in Equation 6.6 of [1].
%
%REFERENCES:
%[1] S.-D. Shih, "The period of a Lotka-Volterra system," Taiwanese Journal
%    of Mathematics, vol. 1, no. 4, pp. 451-470, Dec. 1997.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
val=Phi(s/d).*Phi((E-s)/a);
val(~isfinite(val))=0;%If the endpoints were evaluated.

end

function val=Phi(s)
%%PHI The Phi function of Equation 6.7 of [1].
%
%REFERENCES:
%[1] S.-D. Shih, "The period of a Lotka-Volterra system," Taiwanese Journal
%    of Mathematics, vol. 1, no. 4, pp. 451-470, Dec. 1997.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

z=-exp(-1-s);
val=1./(1+LambW(z,0))-1./(1+LambW(z,-1));
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
