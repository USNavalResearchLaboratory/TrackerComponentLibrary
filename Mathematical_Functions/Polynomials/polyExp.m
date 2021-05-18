function B=polyExp(A,P,direction)
%%POLYEXP Given a polynomial of the form
%         f(x)=A(1)*x^N+A(2)*x^(N-1)+...A(N-1)*x+A(N+1) where A(N+1)=1
%         (if direction=0) or of the form
%         f(x)=A(1)+A(2)*x+A(3)*x^2+...+A(N+1)*x^N where A(1)=1
%         (if direction=1), compute the first N coefficients of a
%         polynomial expansion of f(x)^P, where P is any real number. The
%         polynomial expansion g has the same form as f with 
%         g(x)=B(1)*x^N+B(2)*x^(N-1)+...B(N-1)*x+B(N+1) where B(N+1)=1
%         (if direction=0) or of the form
%         f(x)=B(1)+B(2)*x+B(3)*x^2+...+B(N+1)*x^N where B(1)=1
%         (if direction=1). P=0 is a special case and it results in the
%         natural logarithm being obtained.
%
%INPUTS: A An NX1 or 1XN set of coefficients. Zero-pad the higher-order
%          coefficients to get a higher-order series approximation on the
%          output.
%        P The real power to which the polynomial sum should be raised. P=0
%          is a special case meaning that the logarithm of the series is
%          desired.
% direction An optional parameter specifying which of the two polynomial
%          formats is used. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) A(1)is
%            the coefficient of the hihest order x term.
%          1 A(1) is the coefficient of the lowest order x term.
%
%OUTPUTS: B An NX1 vector of coefficients approximating the series given in
%           A raised to P.
%
%The algorithms is from [2], which corrects the errors in the original
%algorithm in [1] and [3]. In [1], only a minor modification of [3] is made
%so that the natural logarithm is taken for the P=0 case. 
%
%EXAMPLE 1:
%This is Example 1 from [2]. The exact solution is easy to find, because it
%is just squaring a polynomial.
% P=2;
% f=[1/2;3;2;1];
% %By explicitely squaring the polynomial, we can get the exact solution.
% %Convolution is scalar polynomial multiplication:
% f2Exact=conv(f,f);
% %We actually use a zero-padded version of f so as to guarantee having
% %enough coefficients to accurately represent the squared polynomial.
% fPadded=[zeros(3,1);f];
% f2=polyExp(fPadded,P);
% max(abs((f2-f2Exact)./f2Exact))%The error should be 0.
%
%EXAMPLE 2:
%In some instances, particularly with negative, zero or with non-integer P,
%the approximation might only be valid in a limited region about the
%origin. Here, we demonstrate this by plotting the result for P=0, which
%corresponds to taking the natural logarithm of f(x). 
% P=0;
% fPadded=[zeros(20,1);f];%Pad to a higher order.
% fExp=polyExp(fPadded,P);
% 
% x=linspace(-1,1,5000);
% figure(1)
% clf
% hold on
% plot(x,log(polyval(f,x)),'-k','linewidth',4)
% plot(x,polyval(fExp,x),'-r','linewidth',2)
% axis([-1,1,-1,2])
% legend('Exact Logarithm','Series Approximation','location','southwest')
% h1=xlabel('x');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the approximation appears fits well near the origin, but
%diverges around +/-3/4.
%
%EXAMPLE 3:
%This is example 5 from [2]. f is a series expansion of e^x. For this series,
%choosing P=0 gives us log(exp(x))=x.
% P=0;
% f=[1./factorial(10:-1:1),1];
% fP=polyExp(f,P);
% x=linspace(-10,10,500);
% figure(1)
% clf
% hold on
% plot(x,x,'-k','linewidth',4)
% plot(x,polyval(fP,x),'-r','linewidth',2)
% legend('Exact Logarithm','Series Approximation','location','southwest')
% h1=xlabel('x');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
%One can see that the correct linear solution is obtained.
%
%REFERENCES:
%[1] H. E. Fettis, "Algorithm 158 (Algorithm 134 revised): Exponentiation
%    of series," Communications of the ACM, vol. 6, no. 3, pg. 104, Mar.
%    1963.
%[2] J. D. Lawrence, "Certification of Algorithm 158 exponentiation of
%    series," Communications of the ACM, vol. 6, no. 9, pg. 522, Sep. 1963.
%[3] H. E. Fettis, "Algorithm 134: Exponentiation of series,"
%    Communications of the ACM, vol. 5, no. 11, pg. 553, Nov. 1962.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(direction))
    direction=0; 
end

if(direction==0)
    A=flipud(A(:));
end

if(A(1)~=1)
    error('This function only works with series whose constant term is 1.') 
end
N=length(A)-1;
B=zeros(N+1,1);


if(P==0)
    %The constant term is zero when we we are performing a logarithm.
    B(1)=0;
    R=1;
else
    B(1)=1;%The constant term is 1.
    R=P;
end
B(1+1)=R*A(1+1);

for i=2:N
    S=0;
    for k=1:(i-1)
        S=S+(P*(i-k)-k)*B(k+1)*A(i-k+1);
    end
    B(i+1)=R*A(i+1)+(S/i);
end

if(direction==0)
    B=flipud(B);
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
