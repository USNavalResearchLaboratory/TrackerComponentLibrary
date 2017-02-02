function [minPoint,minValue,exitCode]=goldenSectionSearch(g,xSpan,XTol,maxIter)
%%GOLDENSECTIONSEARCH Perform a golden section search to find the minimum
%                     of a continuous, real, function that is unimodal in a
%                     given interval. The golden section search linearly
%                     decreases the search region. The algorithm is
%                     considered more efficient than a ternary search.
%                     Unlike Matlab's built-in fminbnd function, quadratic
%                     interpolation is not used to try to speed up the
%                     convergence rate.
%
%INPUTS: g The handle to a real function to minimize. The function takes a
%          scalar input and returns a real scalar output.
%    xSpan The 2X1 or 1X2 span of values in which the minimum is to be
%          found. xSpan(1)<xSpan(2).
%     XTol The absolute tolerance on the argument of g. When the region of
%          uncertainty of the minimum goes below xTol, then convergence is
%          declared. If omitted or an empty matrix is passed, XTol=1e-9 is
%          used.
%  maxIter The maximum number of iterations to use. If this parameter is
%          omitted or an empty matrix is passed, then maxIter=100; is used.
%          One to 2 function evaluations are performed each iteration.
%
%OUTPUTS: minPoint The value such that g(minVal) is minimized after the
%                  search.
%         minValue The value g(minVal).
%         exitCode A value indicating how the function terminated. Possible
%                  values are:
%                  0 The tolerance XTol  in the argument of g was achieved.
%                  1 The maximum number of iterations was reached.
%
%The golden section search algorithm is implemented as described in
%Appendix C of [1]. The convergence rate is linear.
%
%As an example, consdier finding the minimum of the function
%g(x)=(x-3)*x^3*(x-6)^4.
%This function has multiple minima and maxima. However, if one minimum can
%be bounded, then the golden section search can be used to find it. If
%multiple minima are bounded, then the golden section search can be used,
%but it might not convergence to the best minimum. We will consider the
%the case of three ranges.
% %Example 1
% alphaSpan=[0;3];
% g=@(x)(x-3)*x^3*(x-6)^4;
% [minPoint,minValue]=goldenSectionSearch(g,alphaSpan);
% %Here we find a minPoint of about 1.7354 and a min value of about -2186.1.
% %Next, we consider bounding a second minimum.
% %Example 2
% alphaSpan=[3;7];
% [minPoint1,minValue1]=goldenSectionSearch(g,alphaSpan);
% %Here, we get a minimum point of 6 with a value of 0. Now, if we take a
% %span that covers both minimum:
% %Example 3:
% alphaSpan=[1;7];
% [minPoint2,minValue2]=goldenSectionSearch(g,alphaSpan);
% %We also get a minimum point of 6 even though that is not the absolute
% %minimum over the entire region. This is because the entire region is not
% %unimodal. If we used alphaSpan=[0;7]; then the first minimum (the lower
% one) in that span would have been found.
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 2nd ed. Belmont, MA: Athena
%    Science, 1999.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(maxIter))
    maxIter=100;
end

if(nargin<3||isempty(XTol))
   XTol=1e-9; 
end

tau=(3-sqrt(5))/2;%The golden ratio

alpha=xSpan(1);
alphaBar=xSpan(2);

if(alphaBar<alpha)
   error('alphaSpan(2) is less than alphaSpan(1)');
end

b=alpha+tau*(alphaBar-alpha);
bBar=alphaBar-tau*(alphaBar-alpha);

gAlpha=g(alpha);
gAlphaBar=g(alphaBar);

gb=g(b);
gbBar=g(bBar);

exitCode=1;
for curIter=1:maxIter
    if(gb<gbBar)
        if(gAlpha<=gb)
            alphaBar=b;
            gAlphaBar=gb;
            
            bBar=alphaBar-tau*(alphaBar-alpha);
            gbBar=g(bBar);
        else
            alphaBar=bBar;
            gAlphaBar=gbBar;
            
            %Because of the use of the Golden ratio for tau, only one
            %function evaluation is needed in this instance.
            bBar=b;
            gbBar=gb;
        end

        b=alpha+tau*(alphaBar-alpha);
        gb=g(b);
    elseif(gb>gbBar)
        if(gbBar>=gAlphaBar)
            alpha=bBar;
            gAlpha=gbBar;
            
            b=alpha+tau*(alphaBar-alpha);
            gb=g(b);
        else
            alpha=b;
            gAlpha=gb;
            
            %Because of the use of the golden ratio for tau, only one
            %function evaluation is needed in this instance.
            b=bBar;
            gb=gbBar;
        end
        
        bBar=alphaBar-tau*(alphaBar-alpha);
        gbBar=g(bBar);
    else
        alpha=b;
        alphaBar=bBar;
        
        gAlpha=gb;
        gAlphaBar=gbBar;
        
        
        b=alpha+tau*(alphaBar-alpha);
        gb=g(b);
        bBar=alphaBar-tau*(alphaBar-alpha);
        gbBar=g(bBar);
    end
    
    %If the 
    if(alphaBar-alpha<XTol)
        exitCode=1;
        break;
    end
end

%After the iterations, we have the function evaluated at alpha, b, bBar,
%and alphaBar. We will just return the point/ value that is the smallest.

points=[alpha;b;bBar;alphaBar];
values=[gAlpha;gb;gbBar;gAlphaBar];

[values,idx]=sort(values,'ascend');
minPoint=points(idx(1));
minValue=values(1);
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
