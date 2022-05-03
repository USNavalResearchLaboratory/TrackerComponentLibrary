function y=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1,method)
%%EVALCLENSHAWRECURSERIES Use one of two forms of Clenshaw's method to
%           evaluate a series of the form sum_{k=0}^Nc(k+1)*F(k) where F(k)
%           is a function that is subject to the recurrence relation
%           F(k+1)=alpha(k)*F(k)+beta(k)*F(k-1)
%           Clenshaw's method is numerically stabler than directly using
%           the above recursion and explicitly evaluating the sum.
%           Examples of functions satisfying this type of relation are sine,
%           cosine, Legendre polynomials, and Bessel functions, among
%           others. This function is more efficient if alpha and beta are
%           constants.
%
%INPUTS:  c An (N+1)X1 or 1X(N+1) vector of coefficients for the sum.
%  alphaVal, betaVal These are the coefficients in the recursion on F given
%           above. These can either both be constants, or they can
%           both be function handles that take a single parameter, as
%           implied by the functions alpha(k) and beta(k) shown in the
%           above equation. If function handles, they must return scalar
%           values. If constants, they can be scalars or vectors. if
%           vectors, they must be the same dimensionality as F0 and F1.
%    F0, F1 If method=0, then F0=F(0) and F1=F(1). The first two values of
%           F are needed to start the recursion. On the other hand, if
%           method=1, then F0=F(N) and F1=F(N-1) as the recursion goes in
%           the opposite direction. To evaluate multiple function values at
%           once, F0 and F1 can be matrices. Both must be the same size. If
%           these are matrices and alphaVal and betaVal are not scalars,
%           then alphaVal and betaVal must be the same dimensionality as
%           these.
%    method This optionally selects the type of series to use. Possible
%           values are:
%           0 (The default if omitted or an empty matrix is passed) Use the
%             Clenshaw series going backwards from N such that F0=F(0) and
%             F1=F(1).
%           1 Use the Clenshaw series going forward from 0 such that
%             F0=F(N) and F1=F(N-1). This method is generally only
%             beneficial if F(k) is small when k is large and c(k) is
%             small when k is small.
%
%OUTPUTS: y The value of the sum. This has the same size as F0.
%
%Equations for the two approaches to recursive evaluation of Clenshaw-
%style series are given in Chapter 5.5 of [1]. Note that no code is given
%in Chapter 5.5 and no code given elsewhere in the book was used to
%implement this function.
%
%EXAMPLE 1:
%The recurrence for powers of the cosine function is
% cos(n*theta)=2*cos(theta)*cos((n-1)*theta)-cos((n-2)*theta);
%Suppose we want to evaluate the sum
%sum_{k=0}^N c(k+1)*cos(k*theta)
%This example evaluates the sum for a random set of c and theta using both
%methods 0 and 1 and then compares to explicitly doing out the sum.
%Parameters defining the sum
% N=11;
% c=randn(N+1,1);
% theta=2*pi*rand();
% %Parameters that are the same in methods 0 and 1.
% alphaVal=2*cos(theta);
% betaVal=-1;
% %First use method 0
% method=0;
% F0=1;%=cos(0*theta)
% F1=cos(theta);
% y0=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1,method)
% %Next, use method 1
% method=1;
% F0=cos(N*theta);
% F1=cos((N-1)*theta);
% y1=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1,method)
% %Finally, just directly evaluate the sum.
% k=(0:N).';
% y2=sum(c.*cos(k*theta))
%One will see that all of the sums have the same value, within finite
%precision bounds. Generally, one will want to use method 0.
%
%EXAMPLE 2:
%The recurrence for Bessel functions of the first kind in terms of the
%order is
% besselj(n+1,x)=(2*n/x)besselj(n,x)-besselj(n-1,x);
%Suppose we want to evaluate the sum
%sum_{k=0}^N c(k+1)*besselj(k,x)
%This example evaluates the sum for a random set of c and theta using both
%methods 0 and 1 and then compares to explicitly doing out the sum.
%Parameters defining the sum
% N=11;
% c=randn(N+1,1);
% x=10*randn(1,1);
% %Parameters that are the same in methods 0 and 1.
% alphaVal=@(n)(2*n/x);
% betaVal=@(n)(-1);
% %First use method 0
% method=0;
% F0=besselj(0,x);
% F1=besselj(1,x);
% y0=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1,method)
% %Next, use method 1
% method=1;
% F0=besselj(N,x);
% F1=besselj((N-1),x);
% y1=evalClenshawRecurSeries(c,alphaVal,betaVal,F0,F1,method)
% %Finally, just directly evaluate the sum.
% k=(0:N).';
% y2=sum(c.*besselj(k,x))
%One will see that all of the sums have the same value, within finite
%precision bounds. Generally, one will want to use method 0.
%
%REFERENCES:
%[1] W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling,
%    Numerical Recipes in C, 2nd ed. Cambridge University Press, 1992.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<6||isempty(method))
        method=0; 
    end

    N=length(c)-1;

    sizeF0=size(F0);
    y2=zeros(sizeF0);
    y1=zeros(sizeF0);
    
    if(isnumeric(alphaVal))
        %If the alpha and beta terms are constant.
        switch(method)
            case 0%The standard Clenshaw recursion.
                for k=N:-1:1
                    y=alphaVal.*y1+betaVal.*y2+c(k+1);
                    y2=y1;
                    y1=y;
                end
                y=betaVal.*F0.*y2+F1.*y1+F0.*c(0+1);
            case 1%The Clenshaw recursion going in the opposite direction.
                for k=0:(N-1)
                    y=(y2-alphaVal.*y1-c(k+1))./betaVal;
                    y2=y1;
                    y1=y;
                end
                FN=F0;
                FN1=F1;

                y=c(N+1).*FN-betaVal.*FN1.*y1-FN.*y2;
            otherwise
                error('Unknown method specified.')
        end
    else%If the alpha and beta terms are functions.
         switch(method)
            case 0%The standard Clenshaw recursion.
                for k=N:-1:1
                    y=alphaVal(k).*y1+betaVal(k+1).*y2+c(k+1);
                    y2=y1;
                    y1=y;
                end
                y=betaVal(1).*F0.*y2+F1.*y1+F0.*c(0+1);
            case 1%The Clenshaw recursion going in the opposite direction.
                for k=0:(N-1)
                    y=(y2-alphaVal(k).*y1-c(k+1))./betaVal(k+1);
                    y2=y1;
                    y1=y;
                end
                FN=F0;
                FN1=F1;

                y=c(N+1).*FN-betaVal(N).*FN1.*y1-FN.*y2;
            otherwise
                error('unKnown method specified.')
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
