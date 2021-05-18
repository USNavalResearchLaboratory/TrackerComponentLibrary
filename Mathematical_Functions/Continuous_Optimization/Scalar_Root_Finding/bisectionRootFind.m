function [xMid,fMid,exitCode]=bisectionRootFind(f,xSpan,AbsTol,maxIter)
%%BISECTIONROOTFIND Given two points that are of different signs
%          surrounding the root of a scalar funding, find the root of the
%          function via a bisection search.
%
%INPUTS: f The handle to the function whose root is desired. The function
%          is called as f(x), where x is a scalar value.
%    xSpan A 2X1 or 1X2 vector such that
%          sign(f(xSpan(1)))~=sign(f(xSpan(2)))
%   AbsTol The absolute tolerance on the value of f for converence. The
%          default if omitted or an empty matrix is passed is eps().
%  maxIter The maximum number of bisections to perform. The default if
%          omitted or an empty matrix is passed is 100.
%
%OUTPUTS: xMid The point found that is the value of the approximate zero.
%         fMid The value of the function at xMid.
%     exitCode A value indicating how the function terminated. Possible
%              values are:
%              0 Either AbsTol was satisfied or the bounds on the value
%                xMid went below finite precision constraints.
%              1 The maximum number of iterations occurred.
%
%The bisection algorithm is straightforward: Assuming that only a single
%root is surrounded by xSpan, one evaluates the point at the middle of the
%interval. One then discards the point in xSpan having the same sign as the
%point in the middle of the interval. That gives a new interval and the
%process continues until convergence. The method is discussed in [1].
%
%EXAMPLE:
% f=@(x)(x.^2-4);
% xSpan=[0;5];
% [xMid,fMid,exitCode]=bisectionRootFind(f,xSpan)
%In this example convergence is to the exact root, xMid=2.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Bisection." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Bisection.html
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
if(nargin<4||isempty(maxIter))
    maxIter=100;
end

if(nargin<3||isempty(AbsTol))
    AbsTol=eps(); 
end
    
xLeft=xSpan(1);
xRight=xSpan(2);
    
fLeft=f(xLeft);
fRight=f(xRight);

if(sign(fLeft)==sign(fRight))
    error('bisectionRootFind:NoChange','There is no sign change in the initial set. ');
end

xMid=xRight;
fMid=fRight;
    
exitCode=1;
for curIter=1:maxIter
    xMid=(xLeft+xRight)/2;
    fMid=f(xMid);

    if(abs(fMid)<AbsTol||xMid==xRight||xMid==xLeft)
        exitCode=0;
        break; 
    end

    if(sign(fLeft)~=sign(fMid))
        xRight=xMid;
    else
        fLeft=fMid;
        xLeft=xMid;
    end
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
