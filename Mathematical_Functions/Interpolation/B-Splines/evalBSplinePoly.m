function val=evalBSplinePoly(t,x,k,idx)
%%EVALBSPLINEPOLY Given a set of knots in non-decreasing order (some can be
%                 repeated), evaluate the idx-th order-k B-spline
%                 polynomial. This is B_{idx,k}(x). The first subscript
%                 refers to the knots used in determining the polynomial.
%                 Unlike the function evalBSplinePolys, this function can
%                 handle values of x that are outside of the region
%                 t(idx+-1k)<=x<=t(idx+k).
%
%INPUTS: t A 1XnumKnots or numKnotsX1 vector containing the knots. The
%          knots must be sorted in ascending order. It is possible for the
%          knots to be repeated.
%        x The scalar or matrix set of points at which the polynomials
%          should be evaluated.
%        k Optionally, the order of the B-spline polynomials. If this
%          parameter is omitted, then the default of fix(length(t)/2)+1 is
%          used. This is the highest order that can be used with the given
%          number of knots.
%      idx This input specifies where in terms of the knots the B-spline
%          polynomials are centered. This is an integer value from 1 to
%          numKnots-k.
%
%OUTPUTS: val The value of the b-spline B_{idx,k}(x) evaluated at all of
%             the points in x. This has the same size as x.
%
%This function implements the divided difference formulation of the splines
%from Equation 1 in Chapter IX of [1]. A clearer explanation of the
%notation is given in Chapter 2.4.4 of [2], where Equation 2.4.4.3 is the
%divided difference formula.
%
%EXAMPLE:
%Here, we plot all of the b-spline basis functions from example IX.1 of
%Chapter IX of [1]. The functions only sum to one over the region from 1 to
%6. Thus, the sum is not one over the region plotted before x=1.
% t=[0;1;1;3;4;6;6;6];
% k=3;
% numPoints=500;
% x=linspace(0,6,numPoints);
% 
% b1=evalBSplinePoly(t,x,k,1);
% b2=evalBSplinePoly(t,x,k,2);
% b3=evalBSplinePoly(t,x,k,3);
% b4=evalBSplinePoly(t,x,k,4);
% b5=evalBSplinePoly(t,x,k,5);
% figure(1)
% clf
% hold on
% plot(x,b1,'-k','linewidth',2)
% plot(x,b2,'--r','linewidth',2)
% plot(x,b3,'-.g','linewidth',2)
% plot(x,b4,'-b','linewidth',2)
% plot(x,b5,'--c','linewidth',2)
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%[2] J. Stoer and R. Burlisch, Introduction to Numerical Analysis, 2nd ed.
%    New York: Springer-Verlag, 1991.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=zeros(size(x));
numEls=numel(x);
for curIdx=1:numEls
    xCur=x(curIdx);
    tSel=t(idx:(idx+k));
    f=max(tSel-xCur,0).^(k-1);
    fDeriv=@(ti,nd)prod((k-(1:nd)))*max(ti-xCur,0).^(k-1-nd);

    val(curIdx)=(t(idx+k)-t(idx))*evalDividedDiff(tSel,f,fDeriv);
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
