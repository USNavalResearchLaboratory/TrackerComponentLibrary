function [v,beta]=HouseholderVec(x,forceSign)
%%HOUSEHOLDERVEC Given an mX1 column vector x, compute mX1 column
%          vector v and scalar beta such that P=eye(m,m)-beta*v*v' is
%          orthogonal and P*x=(+/-)norm(x)*e1 for real x, where e1 is a
%          vector with a 1 in the first entry and zeros elsewhere. If
%          forceSign is true, then it will be such that the + sign is
%          always chosen. For complex x,
%          P*x=(+/-)exp(1j*angle(x(1)))*norm(x)*e1, whereby the sign in
%          front can vary. The v vector is chosen such that v(1)=1. If a
%          scalar value is provided, then v=1, beta=0 are returned.
%
%INPUTS: x An mX1 vector, real or complex.
% forceSign If true and x is real, then P*x always equals +norm(x)*e1.
%          Otherwise, it can be + or -. Often, assumptions on recovering
%          matrices from factored form implicitly assume that
%          forceSign=false. The default if omitted or an empty matrix is
%          passed is false.
%
%OUTPUTS: v An mX1 vector such that P=eye(m)-beta*v*v' is an orthogonal
%           matrix and P*x having the only nonzero element be the first
%           one. If a scalar input is given, then v=0 is returned.
%      beta The real scalar beta as mentioned with v. beta= 2/(v'*v). If a
%           scalar input is given, then beta=0 is returned.
%
%For real vectors x, the Householder vector algorithm of Section 5.1.1 of
%[1] is used. For complex vectors, the Householder algorithm of Section
%5.1.13 is used.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(forceSign))
    forceSign=false;
end

m=length(x);

if(isempty(x))
   v=[];
   beta=[];
   return;
end

if(m==1)
   v=1;
   beta=0;
   return;
end

if(isreal(x))%The real case from Chapter 5.1.3 in [1].
    sigma=x(2:m)'*x(2:m);
    v=[1;x(2:m)];
    if(sigma==0&&x(1)>=0)
        if(forceSign)
            beta=0;
        else
            beta=2;
        end

        return;
    elseif(sigma==0&&x(1)<0)
        beta=2;%Sign corrected from that in the book.
        return;
    else
        mu=sqrt(x(1).*x(1)+sigma);
        %The book has <=; we are using < as it does not change things and
        %for some data types, such as Intervals, it is easier to overload <
        %than <=.
        if(x(1)<0)
            v(1)=x(1)-mu;
        else
            v(1)=-sigma/(x(1)+mu);
        end
        beta=2*v(1).*v(1)./(sigma+v(1).*v(1));
        v=v./v(1);
    end
else%The complex case from Chapter 5.1.13
    angVal=x(1)/abs(x(1));
    %Deal with NaNs due to x(1)=0.
    if(~isfinite(angVal))
        angVal=1;
    end

    x2Norm=norm(x,2);
    
    v=zeros(m,1);
    
    v(2:end)=x(2:end);
    v1=x(1)+angVal*x2Norm;
    v2=x(1)-angVal*x2Norm;
    
    %Of v1 and v2, choose the one with the largest norm.
    if(abs(v1)>abs(v2))
        v(1)=v1;
    else
        v(1)=v2;
    end
    
    v=v/v(1);
    beta=2/(v'*v);
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
