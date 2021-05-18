function B=evalBSplinePolys(t,x,k,idx)
%%EVALBSPLINEPOLYS Given a set of at least 2*(k-1) knots in non-decreasing
%                  order (some can be repeated), evaluate all associated
%                  scalar order-k b-spline polynomials that might be
%                  nonzero. These are B_{idx-k+1,k}(x) to B_{idx,k}(x). The
%                  first subscript refers to the knots used in determining
%                  the polynomial and the second subscript refers to the
%                  order of the polynomial. The center index of the knots
%                  is idx and must be k-1<=idx<length(t)-(k-1). These
%                  polynomial evaluations can be used in B-spline
%                  interpolation. This function only gives accurate results
%                  if t(idx)<=x<=t(idx+1). Outside of that range, consider
%                  the function evalBSplinePoly, which evaluates the
%                  polynomials one at a time.
%
%INPUTS: t A 1XnumKnots or numKnotsX1 vector containing the knots. The
%          knots must be sorted in ascending order. It is possible for the
%          knots to be repeated. However, it is required that
%          t(idx)<t(idx+1). If more than the minimum number of knots needed
%          for a given order k are passed, then the extra knots will be
%          ignored. Which knots are "extra" depends on the idx input. These
%          must be real.
%        x The numPointsX1 or 1XnumPoints set of points at which the
%          polynomials should be evaluated. It is assumed that
%          t(idx)<=x<=t(idx+1) for all points. Outside of this range, the
%          results are invalid. These must be real.
%        k Optionally, the order of the B-spline polynomials. If this
%          parameter is omitted, then the default of fix(length(t)/2)+1 are
%          used. This is the highest order that can be used with the given
%          number of knots. This is a positive real value.
%      idx This optional input specifies where in terms of the knots the
%          B-spline polynomials are centered. The default if omitted or an
%          empty matrix is passed is k-1. If t is length 2*(k-1), then idx
%          must be k-1 or be omitted. In general idx is limited to
%          k-1<=idx<length(t)-(k-1).
%
%OUTPUTS: B A kXnumPoints matrix contianing the B-spline polynomials values
%           B_{idx-k+1,k}(x) to B_{idx,k}(x) evaluated at each value in x.
%           These values are real.
%
%B-spline polynomials are discussed in detail in Chapter IX of [1]. The are
%most simply expressed by the recursion starting with B_{idx,0}(x)=1 for
%t(idx)<=x<=t(idx+1) (and is zero otherwise) and then 
%B_{idx,k}(x)=((x-t(idx))/(t(idx+k)-t(idx)))*B_{idx,k-1}(x)+((t(idx+k+1)-x)/(t(idx+k+1)-t(idx+1)))*B_{idx+1,k-1}(x)
%An algorithm to efficiently implement such a recursion is described in
%Chapter X of [1] and is implemented here. This function has been
%implemented such that when faced with knots such that t(idx)=t(idx+1), all
%of the B values are zero, because the non-finite term is replaced with
%zero. Except in such an instance of repeated knots, sum(B)==1.
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

kKnots=fix(length(t)/2)+1;

if(nargin<3||isempty(k))
    k=kKnots; 
end

if(kKnots<k)
   error('There are not enough points to obtain B-spline polynomials of the desired order.')
end

if(nargin<4||isempty(idx))
   idx=k-1; 
end

if(idx<k-1)
    error('The value of iIdx is too small.') 
elseif(idx>length(t)-(k-1))
    error('The value of iIdx is too large.') 
end

%This is doing if(any(x<t(idx)|x>t(idx+1)) However, finite precision errors
%can make things slightly off and create unnecessary warnings, so the
%comparisons with a tolerance help to avoid that.
epsMult=4;
if(any(-(x-t(idx))>epsMult*eps(x)|(x-t(idx+1))>epsMult*eps(x)))
    warning('Points are outside of the valid region implied by iIdx. Function results can be inaccurate.')
end

x=x(:)';
numPoints=length(x);

B=zeros(k,numPoints);

deltaR=zeros(k-1,numPoints);
deltaL=zeros(k-1,numPoints);

B(1,:)=1;
for j=1:(k-1)
    deltaR(j,:)=t(idx+j)-x;
    deltaL(j,:)=x-t(idx+1-j);
    
    recurVal=zeros(1,numPoints);
    for r=1:j
        term=B(r,:)./(deltaR(r,:)+deltaL(j+1-r,:));
        B(r,:)=recurVal+deltaR(r,:).*term;
       
        recurVal=deltaL(j+1-r,:).*term;
    end
    B(j+1,:)=recurVal;
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
