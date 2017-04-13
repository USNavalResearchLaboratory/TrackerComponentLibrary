function [q,r]=polyDivide(u,v)
%%POLYDIVIDE Given two vectors representing 1D polynomials as power series,
%         divide one by the other (u/v) and determine the remainder. The
%         results are such that u==polySum(conv(q,v),r), where the
%         convolution performs polynomial multiplication. The format of
%         the polynomial vectors is the same as that used in Matlab's
%         polyval function: The coefficients are ordered
%         y(z)=u(end)+u(end-1)*z+u(end-2)*z^2+...
%
%INPUTS: u,v Coefficients for the two power series. The number of terms
%            in the different series can be different. Thus u and v do
%            not have the same length. p1(1) is the coefficient of the
%            highest order term. p(end) is the additive constant term. u
%            and v can be row or column vectors.
%
%OUTPUTS: q The polynomial coefficients of u/v without the remainder.
%         r The remainder polynomial for the division of u/v.
%
%This implements Algorithm D in Chapter 4.6.1 of [1].
%
%EXAMPLE 1:
% u=[1;0;1;0;10;10;8;2;8];
% v=[3;0;5;0;9;4;8];
% %Multiply u and v to get v
% z=conv(u,v);
% %We divide z by v to get u back. The remainder is zero.
% [q,r]=polyDivide(z,v)
% z==polySum(conv(q,v),r)%This statement is all true
%
%EXAMPLE 2:
%Here, we try to divide u by v, but it is not divisible at all, so we just
%get back
% u=[1;0;1;0;10;10;8;2;8];
% v=[3;0;5;0;9;4;8];
% [q,r]=polyDivide(u,v)
% u-polySum(conv(q,v),r)%This is all zeros and a value near eps()
%                       %(finite-precision equality).
%
%EXAMPLE 3:
%Here, we have again no remainder.
% a=[1;2];
% b=[4;16];
% c=[1;12;32];
% z=conv(conv(a,b),c);%multiply polynomials a, b, and c
% [q,r]=polyDivide(z,c)
% z==polySum(conv(q,c),r)%This statement is all true
%
%EXAMPLE 4:
%Here, we have a remainder.
% a=[1;2];
% b=[4;16];
% c=[1;12;32];
% z=conv(conv(a,b),c);%Multiply polynomials a, b, and c
% [q,r]=polyDivide(z,c+1)
% z==polySum(conv(q,c+1),r)%This statement is all true
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Seminumerical Algorithms,
%    2nd ed. Reading, MA: Addison-Wesley, 1998, vol. 2.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

u=u(:);
v=v(:);

m=length(u)-1;
n=length(v)-1;

if(m<n)
    q=[];
    r=u;
    return;
end

qDeg=m-n;
q=zeros(m-n+1,1);

for k=(m-n):-1:0
    q(qDeg-k+1)=u(m-(n+k)+1)/v(n-n+1);
    for j=(n+k-1):-1:k 
        u(m-j+1)=u(m-j+1)-q(qDeg-k+1)*v(n-(j-k)+1);
    end
end

r=u(((m-(n-1)):(m-0))+1);

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
