function val=arithmeticGeometricMean(a,b)
%%ARITHMETICGEOMETRICMEAN Compute the arithmetic-geometric mean of two
%                   numbers. The arithmetic-geometric mean is used in a
%                   number of algorithms for computing transcendental
%                   functions such as complete elliptic integrals as
%                   described in [2].
%
%INPUTS: a,b Two real, positive scalar numbers whose arithmetic-geometric
%            mean is desired.
%
%OUTPUTS: val The arithmetic-geometric mean of a and b.
%
%The arithmetic-geometric mean is the iteration
%a_{n+1}=(a_n+b_n)/2; 
%b_{n+1}=sqrt(a_n*b_n);
%Until convergence as given in [1].
%
%EXAMPLE 1:
%As noted at the end of Chapter 2.4 of [2], the quantity
%val=gamma(3/4)^4/arithmeticGeometricMean(1,1/sqrt(2))^2
%equals pi.
%
%EXAMPLE 2:
%As given in Section 1.2 of [2],
%1/arithmeticGeometricMean(1,x) is equal to
%(2/pi)*integral_0^{pi/2}1/sqrt(1-(1-x^2)*sin(theta)^2) dtheta
%Simplifying one finds that
%pi*sqrt(1/(1-x))/(2*arithmeticGeometricMean(1,sqrt(1/(1-x))))
%is equal to ellipke(x) for 0<x<1.
%Thus
% x=0.75;
% val1=pi*sqrt(1/(1-x))/(2*arithmeticGeometricMean(1,sqrt(1/(1-x))))
% val2=ellipke(x)
%One will find that val1=val2.
%
%REFERENCES:
%[1] Abramowitz, M. and Stegun, I. A. (Eds.). "The process of the
%    arithmetic-geometric mean." in Ch. 17.6 in Handbook of Mathematical
%    Functions with Formulas, Graphs, and Mathematical Tables, 9th
%    printing. New York: Dover, pg. 598, 1972.
%[2] J. M. Borwein and B. P. B., Pi and the AGM. New York: John Wiley and
%    Sons, 1987.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Absolute an relative tolerances to obtain convergence to numerical
%precision limits.
AbsTol=eps(max(a,b));
RelTol=eps(1);

a0=a;
b0=b;
while(1)
    a1=(a0+b0)/2; 
    b1=sqrt(a0*b0);
    
    diffMag=abs(a1-b1);
    relMag=min(a1,b1);
    
    if(diffMag<AbsTol||diffMag<RelTol*relMag)
        break;
    end
    
    %Test to see if convergence failed. This can happen if the tolerances
    %are too small.
    if(a0==a1&&b0==b1)
       error('Convergence failed') 
    end
    
    a0=a1;
    b0=b1; 
end
val=a1;
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
